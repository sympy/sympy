# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import abc
import contextlib
import math
import time
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from enum import IntEnum
from random import Random
from sys import float_info
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    NoReturn,
    Optional,
    TypedDict,
    TypeVar,
    Union,
    cast,
)

import attr

from hypothesis.errors import ChoiceTooLarge, Frozen, InvalidArgument, StopTest
from hypothesis.internal.cache import LRUCache
from hypothesis.internal.compat import add_note, floor, int_from_bytes, int_to_bytes
from hypothesis.internal.conjecture.choice import (
    BooleanKWargs,
    BytesKWargs,
    ChoiceKwargsT,
    ChoiceNameT,
    ChoiceT,
    FloatKWargs,
    IntegerKWargs,
    StringKWargs,
    choice_equal,
    choice_from_index,
    choice_key,
    choice_kwargs_equal,
    choice_kwargs_key,
    choice_permitted,
)
from hypothesis.internal.conjecture.floats import float_to_lex, lex_to_float
from hypothesis.internal.conjecture.junkdrawer import (
    IntList,
    gc_cumulative_time,
    uniform,
)
from hypothesis.internal.conjecture.utils import (
    INT_SIZES,
    INT_SIZES_SAMPLER,
    Sampler,
    calc_label_from_name,
    many,
)
from hypothesis.internal.escalation import InterestingOrigin
from hypothesis.internal.floats import (
    SIGNALING_NAN,
    SMALLEST_SUBNORMAL,
    float_to_int,
    int_to_float,
    make_float_clamper,
    next_down,
    next_up,
    sign_aware_lte,
)
from hypothesis.internal.intervalsets import IntervalSet
from hypothesis.reporting import debug_report

if TYPE_CHECKING:
    from typing import TypeAlias

    from typing_extensions import dataclass_transform

    from hypothesis.strategies import SearchStrategy
    from hypothesis.strategies._internal.strategies import Ex
else:
    TypeAlias = object

    def dataclass_transform():
        def wrapper(tp):
            return tp

        return wrapper


TOP_LABEL = calc_label_from_name("top")
TargetObservations = dict[str, Union[int, float]]

T = TypeVar("T")

# index, ir_type, kwargs, forced
MisalignedAt: TypeAlias = tuple[int, ChoiceNameT, ChoiceKwargsT, Optional[ChoiceT]]


class ExtraInformation:
    """A class for holding shared state on a ``ConjectureData`` that should
    be added to the final ``ConjectureResult``."""

    def __repr__(self) -> str:
        return "ExtraInformation({})".format(
            ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items()),
        )

    def has_information(self) -> bool:
        return bool(self.__dict__)


class Status(IntEnum):
    OVERRUN = 0
    INVALID = 1
    VALID = 2
    INTERESTING = 3

    def __repr__(self) -> str:
        return f"Status.{self.name}"


@dataclass_transform()
@attr.s(slots=True, frozen=True)
class StructuralCoverageTag:
    label: int = attr.ib()


STRUCTURAL_COVERAGE_CACHE: dict[int, StructuralCoverageTag] = {}


def structural_coverage(label: int) -> StructuralCoverageTag:
    try:
        return STRUCTURAL_COVERAGE_CACHE[label]
    except KeyError:
        return STRUCTURAL_COVERAGE_CACHE.setdefault(label, StructuralCoverageTag(label))


NASTY_FLOATS = sorted(
    [
        0.0,
        0.5,
        1.1,
        1.5,
        1.9,
        1.0 / 3,
        10e6,
        10e-6,
        1.175494351e-38,
        next_up(0.0),
        float_info.min,
        float_info.max,
        3.402823466e38,
        9007199254740992,
        1 - 10e-6,
        2 + 10e-6,
        1.192092896e-07,
        2.2204460492503131e-016,
    ]
    + [2.0**-n for n in (24, 14, 149, 126)]  # minimum (sub)normals for float16,32
    + [float_info.min / n for n in (2, 10, 1000, 100_000)]  # subnormal in float64
    + [math.inf, math.nan] * 5
    + [SIGNALING_NAN],
    key=float_to_lex,
)
NASTY_FLOATS = list(map(float, NASTY_FLOATS))
NASTY_FLOATS.extend([-x for x in NASTY_FLOATS])

# These caches, especially the kwargs cache, can be quite hot and so we prefer
# LRUCache over LRUReusedCache for performance. We lose scan resistance, but
# that's probably fine here.
FLOAT_INIT_LOGIC_CACHE = LRUCache(4096)
POOLED_KWARGS_CACHE = LRUCache(4096)

COLLECTION_DEFAULT_MAX_SIZE = 10**10  # "arbitrarily large"


class Example:
    """Examples track the hierarchical structure of draws from the byte stream,
    within a single test run.

    Examples are created to mark regions of the byte stream that might be
    useful to the shrinker, such as:
    - The bytes used by a single draw from a strategy.
    - Useful groupings within a strategy, such as individual list elements.
    - Strategy-like helper functions that aren't first-class strategies.
    - Each lowest-level draw of bits or bytes from the byte stream.
    - A single top-level example that spans the entire input.

    Example-tracking allows the shrinker to try "high-level" transformations,
    such as rearranging or deleting the elements of a list, without having
    to understand their exact representation in the byte stream.

    Rather than store each ``Example`` as a rich object, it is actually
    just an index into the ``Examples`` class defined below. This has two
    purposes: Firstly, for most properties of examples we will never need
    to allocate storage at all, because most properties are not used on
    most examples. Secondly, by storing the properties as compact lists
    of integers, we save a considerable amount of space compared to
    Python's normal object size.

    This does have the downside that it increases the amount of allocation
    we do, and slows things down as a result, in some usage patterns because
    we repeatedly allocate the same Example or int objects, but it will
    often dramatically reduce our memory usage, so is worth it.
    """

    __slots__ = ("index", "owner")

    def __init__(self, owner: "Examples", index: int) -> None:
        self.owner = owner
        self.index = index

    def __eq__(self, other: object) -> bool:
        if self is other:
            return True
        if not isinstance(other, Example):
            return NotImplemented
        return (self.owner is other.owner) and (self.index == other.index)

    def __ne__(self, other: object) -> bool:
        if self is other:
            return False
        if not isinstance(other, Example):
            return NotImplemented
        return (self.owner is not other.owner) or (self.index != other.index)

    def __repr__(self) -> str:
        return f"examples[{self.index}]"

    @property
    def label(self) -> int:
        """A label is an opaque value that associates each example with its
        approximate origin, such as a particular strategy class or a particular
        kind of draw."""
        return self.owner.labels[self.owner.label_indices[self.index]]

    @property
    def parent(self) -> Optional[int]:
        """The index of the example that this one is nested directly within."""
        if self.index == 0:
            return None
        return self.owner.parentage[self.index]

    @property
    def ir_start(self) -> int:
        return self.owner.ir_starts[self.index]

    @property
    def ir_end(self) -> int:
        return self.owner.ir_ends[self.index]

    @property
    def depth(self) -> int:
        """Depth of this example in the example tree. The top-level example has a
        depth of 0."""
        return self.owner.depths[self.index]

    @property
    def discarded(self) -> bool:
        """True if this is example's ``stop_example`` call had ``discard`` set to
        ``True``. This means we believe that the shrinker should be able to delete
        this example completely, without affecting the value produced by its enclosing
        strategy. Typically set when a rejection sampler decides to reject a
        generated value and try again."""
        return self.index in self.owner.discarded

    @property
    def ir_length(self) -> int:
        """The number of ir nodes in this example."""
        return self.ir_end - self.ir_start

    @property
    def children(self) -> "list[Example]":
        """The list of all examples with this as a parent, in increasing index
        order."""
        return [self.owner[i] for i in self.owner.children[self.index]]


class ExampleProperty:
    """There are many properties of examples that we calculate by
    essentially rerunning the test case multiple times based on the
    calls which we record in ExampleRecord.

    This class defines a visitor, subclasses of which can be used
    to calculate these properties.
    """

    def __init__(self, examples: "Examples"):
        self.example_stack: "list[int]" = []
        self.examples = examples
        self.example_count = 0
        self.ir_node_count = 0
        self.result: Any = None

    def run(self) -> Any:
        """Rerun the test case with this visitor and return the
        results of ``self.finish()``."""
        self.begin()
        for record in self.examples.trail:
            if record == IR_NODE_RECORD:
                self.ir_node_count += 1
            elif record >= START_EXAMPLE_RECORD:
                self.__push(record - START_EXAMPLE_RECORD)
            else:
                assert record in (
                    STOP_EXAMPLE_DISCARD_RECORD,
                    STOP_EXAMPLE_NO_DISCARD_RECORD,
                )
                self.__pop(discarded=record == STOP_EXAMPLE_DISCARD_RECORD)
        return self.finish()

    def __push(self, label_index: int) -> None:
        i = self.example_count
        assert i < len(self.examples)
        self.start_example(i, label_index=label_index)
        self.example_count += 1
        self.example_stack.append(i)

    def __pop(self, *, discarded: bool) -> None:
        i = self.example_stack.pop()
        self.stop_example(i, discarded=discarded)

    def begin(self) -> None:
        """Called at the beginning of the run to initialise any
        relevant state."""
        self.result = IntList.of_length(len(self.examples))

    def start_example(self, i: int, label_index: int) -> None:
        """Called at the start of each example, with ``i`` the
        index of the example and ``label_index`` the index of
        its label in ``self.examples.labels``."""

    def stop_example(self, i: int, *, discarded: bool) -> None:
        """Called at the end of each example, with ``i`` the
        index of the example and ``discarded`` being ``True`` if ``stop_example``
        was called with ``discard=True``."""

    def finish(self) -> Any:
        return self.result


def calculated_example_property(cls: type[ExampleProperty]) -> Any:
    """Given an ``ExampleProperty`` as above we use this decorator
    to transform it into a lazy property on the ``Examples`` class,
    which has as its value the result of calling ``cls.run()``,
    computed the first time the property is accessed.

    This has the slightly weird result that we are defining nested
    classes which get turned into properties."""
    name = cls.__name__
    cache_name = "__" + name

    def lazy_calculate(self: "Examples") -> Any:
        result = getattr(self, cache_name, None)
        if result is None:
            result = cls(self).run()
            setattr(self, cache_name, result)
        return result

    lazy_calculate.__name__ = cls.__name__
    lazy_calculate.__qualname__ = cls.__qualname__
    return property(lazy_calculate)


STOP_EXAMPLE_DISCARD_RECORD = 1
STOP_EXAMPLE_NO_DISCARD_RECORD = 2
START_EXAMPLE_RECORD = 3

IR_NODE_RECORD = calc_label_from_name("ir draw record")


class ExampleRecord:
    """Records the series of ``start_example``, ``stop_example``, and
    ``draw_bits`` calls so that these may be stored in ``Examples`` and
    replayed when we need to know about the structure of individual
    ``Example`` objects.

    Note that there is significant similarity between this class and
    ``DataObserver``, and the plan is to eventually unify them, but
    they currently have slightly different functions and implementations.
    """

    def __init__(self) -> None:
        self.labels: list[int] = []
        self.__index_of_labels: "dict[int, int] | None" = {}
        self.trail = IntList()
        self.ir_nodes: list[IRNode] = []

    def freeze(self) -> None:
        self.__index_of_labels = None

    def record_ir_draw(self) -> None:
        self.trail.append(IR_NODE_RECORD)

    def start_example(self, label: int) -> None:
        assert self.__index_of_labels is not None
        try:
            i = self.__index_of_labels[label]
        except KeyError:
            i = self.__index_of_labels.setdefault(label, len(self.labels))
            self.labels.append(label)
        self.trail.append(START_EXAMPLE_RECORD + i)

    def stop_example(self, *, discard: bool) -> None:
        if discard:
            self.trail.append(STOP_EXAMPLE_DISCARD_RECORD)
        else:
            self.trail.append(STOP_EXAMPLE_NO_DISCARD_RECORD)


class Examples:
    """A lazy collection of ``Example`` objects, derived from
    the record of recorded behaviour in ``ExampleRecord``.

    Behaves logically as if it were a list of ``Example`` objects,
    but actually mostly exists as a compact store of information
    for them to reference into. All properties on here are best
    understood as the backing storage for ``Example`` and are
    described there.
    """

    def __init__(self, record: ExampleRecord) -> None:
        self.trail = record.trail
        self.labels = record.labels
        self.__length = self.trail.count(
            STOP_EXAMPLE_DISCARD_RECORD
        ) + record.trail.count(STOP_EXAMPLE_NO_DISCARD_RECORD)
        self.__children: "list[Sequence[int]] | None" = None

    class _ir_starts_and_ends(ExampleProperty):
        def begin(self) -> None:
            self.starts = IntList.of_length(len(self.examples))
            self.ends = IntList.of_length(len(self.examples))

        def start_example(self, i: int, label_index: int) -> None:
            self.starts[i] = self.ir_node_count

        def stop_example(self, i: int, *, discarded: bool) -> None:
            self.ends[i] = self.ir_node_count

        def finish(self) -> tuple[IntList, IntList]:
            return (self.starts, self.ends)

    ir_starts_and_ends: "tuple[IntList, IntList]" = calculated_example_property(
        _ir_starts_and_ends
    )

    @property
    def ir_starts(self) -> IntList:
        return self.ir_starts_and_ends[0]

    @property
    def ir_ends(self) -> IntList:
        return self.ir_starts_and_ends[1]

    class _discarded(ExampleProperty):
        def begin(self) -> None:
            self.result: set[int] = set()

        def finish(self) -> frozenset[int]:
            return frozenset(self.result)

        def stop_example(self, i: int, *, discarded: bool) -> None:
            if discarded:
                self.result.add(i)

    discarded: frozenset[int] = calculated_example_property(_discarded)

    class _parentage(ExampleProperty):
        def stop_example(self, i: int, *, discarded: bool) -> None:
            if i > 0:
                self.result[i] = self.example_stack[-1]

    parentage: IntList = calculated_example_property(_parentage)

    class _depths(ExampleProperty):
        def begin(self) -> None:
            self.result = IntList.of_length(len(self.examples))

        def start_example(self, i: int, label_index: int) -> None:
            self.result[i] = len(self.example_stack)

    depths: IntList = calculated_example_property(_depths)

    class _label_indices(ExampleProperty):
        def start_example(self, i: int, label_index: int) -> None:
            self.result[i] = label_index

    label_indices: IntList = calculated_example_property(_label_indices)

    class _mutator_groups(ExampleProperty):
        def begin(self) -> None:
            self.groups: "dict[int, set[tuple[int, int]]]" = defaultdict(set)

        def start_example(self, i: int, label_index: int) -> None:
            # TODO should we discard start == end cases? occurs for eg st.data()
            # which is conditionally or never drawn from. arguably swapping
            # nodes with the empty list is a useful mutation enabled by start == end?
            key = (self.examples[i].ir_start, self.examples[i].ir_end)
            self.groups[label_index].add(key)

        def finish(self) -> Iterable[set[tuple[int, int]]]:
            # Discard groups with only one example, since the mutator can't
            # do anything useful with them.
            return [g for g in self.groups.values() if len(g) >= 2]

    mutator_groups: list[set[tuple[int, int]]] = calculated_example_property(
        _mutator_groups
    )

    @property
    def children(self) -> list[Sequence[int]]:
        if self.__children is None:
            children = [IntList() for _ in range(len(self))]
            for i, p in enumerate(self.parentage):
                if i > 0:
                    children[p].append(i)
            # Replace empty children lists with a tuple to reduce
            # memory usage.
            for i, c in enumerate(children):
                if not c:
                    children[i] = ()  # type: ignore
            self.__children = children  # type: ignore
        return self.__children  # type: ignore

    def __len__(self) -> int:
        return self.__length

    def __getitem__(self, i: int) -> Example:
        assert isinstance(i, int)
        n = len(self)
        if i < -n or i >= n:
            raise IndexError(f"Index {i} out of range [-{n}, {n})")
        if i < 0:
            i += n
        return Example(self, i)

    # not strictly necessary as we have len/getitem, but required for mypy.
    # https://github.com/python/mypy/issues/9737
    def __iter__(self) -> Iterator[Example]:
        for i in range(len(self)):
            yield self[i]


class _Overrun:
    status: Status = Status.OVERRUN

    def __repr__(self) -> str:
        return "Overrun"


Overrun = _Overrun()

global_test_counter = 0


MAX_DEPTH = 100


class DataObserver:
    """Observer class for recording the behaviour of a
    ConjectureData object, primarily used for tracking
    the behaviour in the tree cache."""

    def conclude_test(
        self,
        status: Status,
        interesting_origin: Optional[InterestingOrigin],
    ) -> None:
        """Called when ``conclude_test`` is called on the
        observed ``ConjectureData``, with the same arguments.

        Note that this is called after ``freeze`` has completed.
        """

    def kill_branch(self) -> None:
        """Mark this part of the tree as not worth re-exploring."""

    def draw_integer(
        self, value: int, *, kwargs: IntegerKWargs, was_forced: bool
    ) -> None:
        pass

    def draw_float(
        self, value: float, *, kwargs: FloatKWargs, was_forced: bool
    ) -> None:
        pass

    def draw_string(
        self, value: str, *, kwargs: StringKWargs, was_forced: bool
    ) -> None:
        pass

    def draw_bytes(
        self, value: bytes, *, kwargs: BytesKWargs, was_forced: bool
    ) -> None:
        pass

    def draw_boolean(
        self, value: bool, *, kwargs: BooleanKWargs, was_forced: bool
    ) -> None:
        pass


@attr.s(slots=True, repr=False, eq=False)
class IRNode:
    ir_type: ChoiceNameT = attr.ib()
    value: ChoiceT = attr.ib()
    kwargs: ChoiceKwargsT = attr.ib()
    was_forced: bool = attr.ib()
    index: Optional[int] = attr.ib(default=None)

    def copy(
        self,
        *,
        with_value: Optional[ChoiceT] = None,
        with_kwargs: Optional[ChoiceKwargsT] = None,
    ) -> "IRNode":
        # we may want to allow this combination in the future, but for now it's
        # a footgun.
        if self.was_forced:
            assert with_value is None, "modifying a forced node doesn't make sense"
        # explicitly not copying index. node indices are only assigned via
        # ExampleRecord. This prevents footguns with relying on stale indices
        # after copying.
        return IRNode(
            ir_type=self.ir_type,
            value=self.value if with_value is None else with_value,
            kwargs=self.kwargs if with_kwargs is None else with_kwargs,
            was_forced=self.was_forced,
        )

    @property
    def trivial(self) -> bool:
        """
        A node is trivial if it cannot be simplified any further. This does not
        mean that modifying a trivial node can't produce simpler test cases when
        viewing the tree as a whole. Just that when viewing this node in
        isolation, this is the simplest the node can get.
        """
        if self.was_forced:
            return True

        if self.ir_type != "float":
            zero_value = choice_from_index(0, self.ir_type, self.kwargs)
            return choice_equal(self.value, zero_value)
        else:
            kwargs = cast(FloatKWargs, self.kwargs)
            min_value = kwargs["min_value"]
            max_value = kwargs["max_value"]
            shrink_towards = 0.0

            if min_value == -math.inf and max_value == math.inf:
                return choice_equal(self.value, shrink_towards)

            if (
                not math.isinf(min_value)
                and not math.isinf(max_value)
                and math.ceil(min_value) <= math.floor(max_value)
            ):
                # the interval contains an integer. the simplest integer is the
                # one closest to shrink_towards
                shrink_towards = max(math.ceil(min_value), shrink_towards)
                shrink_towards = min(math.floor(max_value), shrink_towards)
                return choice_equal(self.value, float(shrink_towards))

            # the real answer here is "the value in [min_value, max_value] with
            # the lowest denominator when represented as a fraction".
            # It would be good to compute this correctly in the future, but it's
            # also not incorrect to be conservative here.
            return False

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, IRNode):
            return NotImplemented

        return (
            self.ir_type == other.ir_type
            and choice_equal(self.value, other.value)
            and choice_kwargs_equal(self.ir_type, self.kwargs, other.kwargs)
            and self.was_forced == other.was_forced
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.ir_type,
                choice_key(self.value),
                choice_kwargs_key(self.ir_type, self.kwargs),
                self.was_forced,
            )
        )

    def __repr__(self) -> str:
        # repr to avoid "BytesWarning: str() on a bytes instance" for bytes nodes
        forced_marker = " [forced]" if self.was_forced else ""
        return f"{self.ir_type} {self.value!r}{forced_marker} {self.kwargs!r}"


@attr.s(slots=True)
class NodeTemplate:
    type: Literal["simplest"] = attr.ib()
    count: Optional[int] = attr.ib()

    def __attrs_post_init__(self) -> None:
        if self.count is not None:
            assert self.count > 0


def ir_size(ir: Iterable[ChoiceT]) -> int:
    from hypothesis.database import ir_to_bytes

    return len(ir_to_bytes(ir))


@dataclass_transform()
@attr.s(slots=True)
class ConjectureResult:
    """Result class storing the parts of ConjectureData that we
    will care about after the original ConjectureData has outlived its
    usefulness."""

    status: Status = attr.ib()
    interesting_origin: Optional[InterestingOrigin] = attr.ib()
    buffer: bytes = attr.ib()
    ir_nodes: tuple[IRNode, ...] = attr.ib(eq=False, repr=False)
    output: str = attr.ib()
    extra_information: Optional[ExtraInformation] = attr.ib()
    has_discards: bool = attr.ib()
    target_observations: TargetObservations = attr.ib()
    tags: frozenset[StructuralCoverageTag] = attr.ib()
    examples: Examples = attr.ib(repr=False, eq=False)
    arg_slices: set[tuple[int, int]] = attr.ib(repr=False)
    slice_comments: dict[tuple[int, int], str] = attr.ib(repr=False)
    misaligned_at: Optional[MisalignedAt] = attr.ib(repr=False)

    index: int = attr.ib(init=False)

    def __attrs_post_init__(self) -> None:
        self.index = len(self.buffer)

    def as_result(self) -> "ConjectureResult":
        return self

    @property
    def choices(self) -> tuple[ChoiceT, ...]:
        return tuple(node.value for node in self.ir_nodes)


# Masks for masking off the first byte of an n-bit buffer.
# The appropriate mask is stored at position n % 8.
BYTE_MASKS = [(1 << n) - 1 for n in range(8)]
BYTE_MASKS[0] = 255

_Lifetime: TypeAlias = Literal["test_case", "test_function"]


class _BackendInfoMsg(TypedDict):
    type: str
    title: str
    content: Union[str, dict[str, Any]]


class PrimitiveProvider(abc.ABC):
    # This is the low-level interface which would also be implemented
    # by e.g. CrossHair, by an Atheris-hypothesis integration, etc.
    # We'd then build the structured tree handling, database and replay
    # support, etc. on top of this - so all backends get those for free.
    #
    # See https://github.com/HypothesisWorks/hypothesis/issues/3086

    # How long a provider instance is used for. One of test_function or
    # test_case. Defaults to test_function.
    #
    # If test_function, a single provider instance will be instantiated and used
    # for the entirety of each test function. I.e., roughly one provider per
    # @given annotation. This can be useful if you need to track state over many
    # executions to a test function.
    #
    # This lifetime will cause None to be passed for the ConjectureData object
    # in PrimitiveProvider.__init__, because that object is instantiated per
    # test case.
    #
    # If test_case, a new provider instance will be instantiated and used each
    # time hypothesis tries to generate a new input to the test function. This
    # lifetime can access the passed ConjectureData object.
    #
    # Non-hypothesis providers probably want to set a lifetime of test_function.
    lifetime: _Lifetime = "test_function"

    # Solver-based backends such as hypothesis-crosshair use symbolic values
    # which record operations performed on them in order to discover new paths.
    # If avoid_realization is set to True, hypothesis will avoid interacting with
    # symbolic choices returned by the provider in any way that would force the
    # solver to narrow the range of possible values for that symbolic.
    #
    # Setting this to True disables some hypothesis features, such as
    # DataTree-based deduplication, and some internal optimizations, such as
    # caching kwargs. Only enable this if it is necessary for your backend.
    avoid_realization = False

    def __init__(self, conjecturedata: Optional["ConjectureData"], /) -> None:
        self._cd = conjecturedata

    def per_test_case_context_manager(self):
        return contextlib.nullcontext()

    def realize(self, value: T) -> T:
        """
        Called whenever hypothesis requires a concrete (non-symbolic) value from
        a potentially symbolic value. Hypothesis will not check that `value` is
        symbolic before calling `realize`, so you should handle the case where
        `value` is non-symbolic.

        The returned value should be non-symbolic.  If you cannot provide a value,
        raise hypothesis.errors.BackendCannotProceed("discard_test_case")
        """
        return value

    def observe_test_case(self) -> dict[str, Any]:
        """Called at the end of the test case when observability mode is active.

        The return value should be a non-symbolic json-encodable dictionary,
        and will be included as `observation["metadata"]["backend"]`.
        """
        return {}

    def observe_information_messages(
        self, *, lifetime: _Lifetime
    ) -> Iterable[_BackendInfoMsg]:
        """Called at the end of each test case and again at end of the test function.

        Return an iterable of `{type: info/alert/error, title: str, content: str|dict}`
        dictionaries to be delivered as individual information messages.
        (Hypothesis adds the `run_start` timestamp and `property` name for you.)
        """
        assert lifetime in ("test_case", "test_function")
        yield from []

    @abc.abstractmethod
    def draw_boolean(
        self,
        p: float = 0.5,
        *,
        forced: Optional[bool] = None,
        fake_forced: bool = False,
    ) -> bool:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_integer(
        self,
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        # weights are for choosing an element index from a bounded range
        weights: Optional[dict[int, float]] = None,
        shrink_towards: int = 0,
        forced: Optional[int] = None,
        fake_forced: bool = False,
    ) -> int:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
        # TODO: consider supporting these float widths at the IR level in the
        # future.
        # width: Literal[16, 32, 64] = 64,
        # exclude_min and exclude_max handled higher up,
        forced: Optional[float] = None,
        fake_forced: bool = False,
    ) -> float:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        forced: Optional[str] = None,
        fake_forced: bool = False,
    ) -> str:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        *,
        forced: Optional[bytes] = None,
        fake_forced: bool = False,
    ) -> bytes:
        raise NotImplementedError

    def span_start(self, label: int, /) -> None:  # noqa: B027  # non-abstract noop
        """Marks the beginning of a semantically meaningful span.

        Providers can optionally track this data to learn which sub-sequences
        of draws correspond to a higher-level object, recovering the parse tree.
        `label` is an opaque integer, which will be shared by all spans drawn
        from a particular strategy.

        This method is called from ConjectureData.start_example().
        """

    def span_end(self, discard: bool, /) -> None:  # noqa: B027  # non-abstract noop
        """Marks the end of a semantically meaningful span.

        `discard` is True when the draw was filtered out or otherwise marked as
        unlikely to contribute to the input data as seen by the user's test.
        Note however that side effects can make this determination unsound.

        This method is called from ConjectureData.stop_example().
        """


class HypothesisProvider(PrimitiveProvider):
    lifetime = "test_case"

    def __init__(self, conjecturedata: Optional["ConjectureData"], /):
        super().__init__(conjecturedata)

    def draw_boolean(
        self,
        p: float = 0.5,
        *,
        forced: Optional[bool] = None,
        fake_forced: bool = False,
    ) -> bool:
        """Return True with probability p (assuming a uniform generator),
        shrinking towards False. If ``forced`` is set to a non-None value, this
        will always return that value but will write choices appropriate to having
        drawn that value randomly."""
        # Note that this could also be implemented in terms of draw_integer().

        assert self._cd is not None
        # NB this function is vastly more complicated than it may seem reasonable
        # for it to be. This is because it is used in a lot of places and it's
        # important for it to shrink well, so it's worth the engineering effort.

        if p <= 0 or p >= 1:
            bits = 1
        else:
            # When there is a meaningful draw, in order to shrink well we will
            # set things up so that 0 and 1 always correspond to False and True
            # respectively. This means we want enough bits available that in a
            # draw we will always have at least one truthy value and one falsey
            # value.
            bits = math.ceil(-math.log(min(p, 1 - p), 2))
        # In order to avoid stupidly large draws where the probability is
        # effectively zero or one, we treat probabilities of under 2^-64 to be
        # effectively zero.
        if bits > 64:
            # There isn't enough precision near one for this to occur for values
            # far from 0.
            p = 0.0
            bits = 1

        size = 2**bits

        while True:
            # The logic here is a bit complicated and special cased to make it
            # play better with the shrinker.

            # We imagine partitioning the real interval [0, 1] into 2**n equal parts
            # and looking at each part and whether its interior is wholly <= p
            # or wholly >= p. At most one part can be neither.

            # We then pick a random part. If it's wholly on one side or the other
            # of p then we use that as the answer. If p is contained in the
            # interval then we start again with a new probability that is given
            # by the fraction of that interval that was <= our previous p.

            # We then take advantage of the fact that we have control of the
            # labelling to make this shrink better, using the following tricks:

            # If p is <= 0 or >= 1 the result of this coin is certain. We make sure
            # to write a byte to the data stream anyway so that these don't cause
            # difficulties when shrinking.
            if p <= 0:
                self._cd.draw_bits(1, forced=0)
                result = False
            elif p >= 1:
                self._cd.draw_bits(1, forced=1)
                result = True
            else:
                falsey = floor(size * (1 - p))
                truthy = floor(size * p)
                remainder = size * p - truthy

                if falsey + truthy == size:
                    partial = False
                else:
                    partial = True

                i = self._cd.draw_bits(
                    bits,
                    forced=None if forced is None else int(forced),
                    fake_forced=fake_forced,
                )

                # We always choose the region that causes us to repeat the loop as
                # the maximum value, so that shrinking the drawn bits never causes
                # us to need to draw more self._cd.
                if partial and i == size - 1:
                    p = remainder
                    continue
                if falsey == 0:
                    # Every other partition is truthy, so the result is true
                    result = True
                elif truthy == 0:
                    # Every other partition is falsey, so the result is false
                    result = False
                elif i <= 1:
                    # We special case so that zero is always false and 1 is always
                    # true which makes shrinking easier because we can always
                    # replace a truthy block with 1. This has the slightly weird
                    # property that shrinking from 2 to 1 can cause the result to
                    # grow, but the shrinker always tries 0 and 1 first anyway, so
                    # this will usually be fine.
                    result = bool(i)
                else:
                    # Originally everything in the region 0 <= i < falsey was false
                    # and everything above was true. We swapped one truthy element
                    # into this region, so the region becomes 0 <= i <= falsey
                    # except for i = 1. We know i > 1 here, so the test for truth
                    # becomes i > falsey.
                    result = i > falsey

            break
        return result

    def draw_integer(
        self,
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        weights: Optional[dict[int, float]] = None,
        shrink_towards: int = 0,
        forced: Optional[int] = None,
        fake_forced: bool = False,
    ) -> int:
        assert self._cd is not None

        if min_value is not None:
            shrink_towards = max(min_value, shrink_towards)
        if max_value is not None:
            shrink_towards = min(max_value, shrink_towards)

        # This is easy to build on top of our existing conjecture utils,
        # and it's easy to build sampled_from and weighted_coin on this.
        if weights is not None:
            assert min_value is not None
            assert max_value is not None

            # format of weights is a mapping of ints to p, where sum(p) < 1.
            # The remaining probability mass is uniformly distributed over
            # *all* ints (not just the unmapped ones; this is somewhat undesirable,
            # but simplifies things).
            #
            # We assert that sum(p) is strictly less than 1 because it simplifies
            # handling forced values when we can force into the unmapped probability
            # mass. We should eventually remove this restriction.
            sampler = Sampler(
                [1 - sum(weights.values()), *weights.values()], observe=False
            )
            # if we're forcing, it's easiest to force into the unmapped probability
            # mass and then force the drawn value after.
            idx = sampler.sample(
                self._cd, forced=None if forced is None else 0, fake_forced=fake_forced
            )

            return self._draw_bounded_integer(
                min_value,
                max_value,
                # implicit reliance on dicts being sorted for determinism
                forced=forced if idx == 0 else list(weights)[idx - 1],
                center=shrink_towards,
                fake_forced=fake_forced,
            )

        if min_value is None and max_value is None:
            return self._draw_unbounded_integer(forced=forced, fake_forced=fake_forced)

        if min_value is None:
            assert max_value is not None  # make mypy happy
            probe = max_value + 1
            while max_value < probe:
                probe = shrink_towards + self._draw_unbounded_integer(
                    forced=None if forced is None else forced - shrink_towards,
                    fake_forced=fake_forced,
                )
            return probe

        if max_value is None:
            assert min_value is not None
            probe = min_value - 1
            while probe < min_value:
                probe = shrink_towards + self._draw_unbounded_integer(
                    forced=None if forced is None else forced - shrink_towards,
                    fake_forced=fake_forced,
                )
            return probe

        return self._draw_bounded_integer(
            min_value,
            max_value,
            center=shrink_towards,
            forced=forced,
            fake_forced=fake_forced,
        )

    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
        # TODO: consider supporting these float widths at the IR level in the
        # future.
        # width: Literal[16, 32, 64] = 64,
        # exclude_min and exclude_max handled higher up,
        forced: Optional[float] = None,
        fake_forced: bool = False,
    ) -> float:
        (
            sampler,
            forced_sign_bit,
            clamper,
            nasty_floats,
        ) = self._draw_float_init_logic(
            min_value=min_value,
            max_value=max_value,
            allow_nan=allow_nan,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
        )

        assert self._cd is not None

        while True:
            # If `forced in nasty_floats`, then `forced` was *probably*
            # generated by drawing a nonzero index from the sampler. However, we
            # have no obligation to generate it that way when forcing. In particular,
            # i == 0 is able to produce all possible floats, and the forcing
            # logic is simpler if we assume this choice.
            forced_i = None if forced is None else 0
            i = (
                sampler.sample(self._cd, forced=forced_i, fake_forced=fake_forced)
                if sampler
                else 0
            )
            if i == 0:
                result = self._draw_float(
                    forced_sign_bit=forced_sign_bit,
                    forced=forced,
                    fake_forced=fake_forced,
                )
                if allow_nan and math.isnan(result):
                    clamped = result
                else:
                    clamped = clamper(result)
                if clamped != result and not (math.isnan(result) and allow_nan):
                    self._draw_float(forced=clamped, fake_forced=fake_forced)
                    result = clamped
            else:
                result = nasty_floats[i - 1]
                # nan values generated via int_to_float break list membership:
                #
                #  >>> n = 18444492273895866368
                # >>> assert math.isnan(int_to_float(n))
                # >>> assert int_to_float(n) not in [int_to_float(n)]
                #
                # because int_to_float nans are not equal in the sense of either
                # `a == b` or `a is b`.
                #
                # This can lead to flaky errors when collections require unique
                # floats. I think what is happening is that in some places we
                # provide math.nan, and in others we provide
                # int_to_float(float_to_int(math.nan)), and which one gets used
                # is not deterministic across test iterations.
                #
                # As a (temporary?) fix, we'll *always* generate nan values which
                # are not equal in the identity sense.
                #
                # see also https://github.com/HypothesisWorks/hypothesis/issues/3926.
                if math.isnan(result):
                    result = int_to_float(float_to_int(result))

                self._draw_float(forced=result, fake_forced=fake_forced)

            return result

    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        forced: Optional[str] = None,
        fake_forced: bool = False,
    ) -> str:
        assert self._cd is not None

        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )

        chars = []
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            forced=None if forced is None else len(forced),
            fake_forced=fake_forced,
            observe=False,
        )
        while elements.more():
            forced_i: Optional[int] = None
            if forced is not None:
                c = forced[elements.count - 1]
                forced_i = intervals.index_from_char_in_shrink_order(c)

            if len(intervals) > 256:
                if self.draw_boolean(
                    0.2,
                    forced=None if forced_i is None else forced_i > 255,
                    fake_forced=fake_forced,
                ):
                    i = self._draw_bounded_integer(
                        256,
                        len(intervals) - 1,
                        forced=forced_i,
                        fake_forced=fake_forced,
                    )
                else:
                    i = self._draw_bounded_integer(
                        0, 255, forced=forced_i, fake_forced=fake_forced
                    )
            else:
                i = self._draw_bounded_integer(
                    0, len(intervals) - 1, forced=forced_i, fake_forced=fake_forced
                )

            chars.append(intervals.char_in_shrink_order(i))

        return "".join(chars)

    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        *,
        forced: Optional[bytes] = None,
        fake_forced: bool = False,
    ) -> bytes:
        assert self._cd is not None

        buf = bytearray()
        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            forced=None if forced is None else len(forced),
            fake_forced=fake_forced,
            observe=False,
        )
        while elements.more():
            forced_i: Optional[int] = None
            if forced is not None:
                # implicit conversion from bytes to int by indexing here
                forced_i = forced[elements.count - 1]

            buf += self._cd.draw_bits(
                8, forced=forced_i, fake_forced=fake_forced
            ).to_bytes(1, "big")

        return bytes(buf)

    def _draw_float(
        self,
        forced_sign_bit: Optional[int] = None,
        *,
        forced: Optional[float] = None,
        fake_forced: bool = False,
    ) -> float:
        """
        Helper for draw_float which draws a random 64-bit float.
        """
        assert self._cd is not None

        if forced is not None:
            # sign_aware_lte(forced, -0.0) does not correctly handle the
            # math.nan case here.
            forced_sign_bit = math.copysign(1, forced) == -1
        is_negative = self._cd.draw_bits(
            1, forced=forced_sign_bit, fake_forced=fake_forced
        )
        f = lex_to_float(
            self._cd.draw_bits(
                64,
                forced=None if forced is None else float_to_lex(abs(forced)),
                fake_forced=fake_forced,
            )
        )
        return -f if is_negative else f

    def _draw_unbounded_integer(
        self, *, forced: Optional[int] = None, fake_forced: bool = False
    ) -> int:
        assert self._cd is not None
        forced_i = None
        if forced is not None:
            # Using any bucket large enough to contain this integer would be a
            # valid way to force it. This is because an n bit integer could have
            # been drawn from a bucket of size n, or from any bucket of size
            # m > n.
            # We'll always choose the smallest eligible bucket here.

            # We need an extra bit to handle forced signed integers. INT_SIZES
            # is interpreted as unsigned sizes.
            bit_size = forced.bit_length() + 1
            size = min(size for size in INT_SIZES if bit_size <= size)
            forced_i = INT_SIZES.index(size)

        size = INT_SIZES[
            INT_SIZES_SAMPLER.sample(self._cd, forced=forced_i, fake_forced=fake_forced)
        ]

        forced_r = None
        if forced is not None:
            forced_r = forced
            forced_r <<= 1
            if forced < 0:
                forced_r = -forced_r
                forced_r |= 1

        r = self._cd.draw_bits(size, forced=forced_r, fake_forced=fake_forced)
        sign = r & 1
        r >>= 1
        if sign:
            r = -r
        return r

    def _draw_bounded_integer(
        self,
        lower: int,
        upper: int,
        *,
        center: Optional[int] = None,
        forced: Optional[int] = None,
        fake_forced: bool = False,
        _vary_effective_size: bool = True,
    ) -> int:
        assert lower <= upper
        assert forced is None or lower <= forced <= upper
        assert self._cd is not None
        if lower == upper:
            # Write a value even when this is trivial so that when a bound depends
            # on other values we don't suddenly disappear when the gap shrinks to
            # zero - if that happens then often the data stream becomes misaligned
            # and we fail to shrink in cases where we really should be able to.
            self._cd.draw_bits(1, forced=0)
            return int(lower)

        if center is None:
            center = lower
        center = min(max(center, lower), upper)

        if center == upper:
            above = False
        elif center == lower:
            above = True
        else:
            force_above = None if forced is None else forced < center
            above = not self._cd.draw_bits(
                1, forced=force_above, fake_forced=fake_forced
            )

        if above:
            gap = upper - center
        else:
            gap = center - lower

        assert gap > 0

        bits = gap.bit_length()
        probe = gap + 1

        if (
            bits > 24
            and _vary_effective_size
            and self.draw_boolean(
                7 / 8, forced=None if forced is None else False, fake_forced=fake_forced
            )
        ):
            # For large ranges, we combine the uniform random distribution from draw_bits
            # with a weighting scheme with moderate chance.  Cutoff at 2 ** 24 so that our
            # choice of unicode characters is uniform but the 32bit distribution is not.
            idx = INT_SIZES_SAMPLER.sample(self._cd)
            force_bits = min(bits, INT_SIZES[idx])
            forced = self._draw_bounded_integer(
                lower=center if above else max(lower, center - 2**force_bits - 1),
                upper=center if not above else min(upper, center + 2**force_bits - 1),
                _vary_effective_size=False,
            )

            assert lower <= forced <= upper

        while probe > gap:
            probe = self._cd.draw_bits(
                bits,
                forced=None if forced is None else abs(forced - center),
                fake_forced=fake_forced,
            )

        if above:
            result = center + probe
        else:
            result = center - probe

        assert lower <= result <= upper
        assert forced is None or result == forced, (result, forced, center, above)
        return result

    @classmethod
    def _draw_float_init_logic(
        cls,
        *,
        min_value: float,
        max_value: float,
        allow_nan: bool,
        smallest_nonzero_magnitude: float,
    ) -> tuple[
        Optional[Sampler],
        Optional[Literal[0, 1]],
        Callable[[float], float],
        list[float],
    ]:
        """
        Caches initialization logic for draw_float, as an alternative to
        computing this for *every* float draw.
        """
        # float_to_int allows us to distinguish between e.g. -0.0 and 0.0,
        # even in light of hash(-0.0) == hash(0.0) and -0.0 == 0.0.
        key = (
            float_to_int(min_value),
            float_to_int(max_value),
            allow_nan,
            float_to_int(smallest_nonzero_magnitude),
        )
        if key in FLOAT_INIT_LOGIC_CACHE:
            return FLOAT_INIT_LOGIC_CACHE[key]

        result = cls._compute_draw_float_init_logic(
            min_value=min_value,
            max_value=max_value,
            allow_nan=allow_nan,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
        )
        FLOAT_INIT_LOGIC_CACHE[key] = result
        return result

    @staticmethod
    def _compute_draw_float_init_logic(
        *,
        min_value: float,
        max_value: float,
        allow_nan: bool,
        smallest_nonzero_magnitude: float,
    ) -> tuple[
        Optional[Sampler],
        Optional[Literal[0, 1]],
        Callable[[float], float],
        list[float],
    ]:
        if smallest_nonzero_magnitude == 0.0:  # pragma: no cover
            raise FloatingPointError(
                "Got allow_subnormal=True, but we can't represent subnormal floats "
                "right now, in violation of the IEEE-754 floating-point "
                "specification.  This is usually because something was compiled with "
                "-ffast-math or a similar option, which sets global processor state.  "
                "See https://simonbyrne.github.io/notes/fastmath/ for a more detailed "
                "writeup - and good luck!"
            )

        def permitted(f: float) -> bool:
            if math.isnan(f):
                return allow_nan
            if 0 < abs(f) < smallest_nonzero_magnitude:
                return False
            return sign_aware_lte(min_value, f) and sign_aware_lte(f, max_value)

        boundary_values = [
            min_value,
            next_up(min_value),
            min_value + 1,
            max_value - 1,
            next_down(max_value),
            max_value,
        ]
        nasty_floats = [f for f in NASTY_FLOATS + boundary_values if permitted(f)]
        weights = [0.2 * len(nasty_floats)] + [0.8] * len(nasty_floats)
        sampler = Sampler(weights, observe=False) if nasty_floats else None

        forced_sign_bit: Optional[Literal[0, 1]] = None
        if sign_aware_lte(min_value, -0.0) != sign_aware_lte(0.0, max_value):
            forced_sign_bit = 1 if sign_aware_lte(min_value, -0.0) else 0

        clamper = make_float_clamper(
            min_value,
            max_value,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
            allow_nan=allow_nan,
        )
        return (sampler, forced_sign_bit, clamper, nasty_floats)


# The set of available `PrimitiveProvider`s, by name.  Other libraries, such as
# crosshair, can implement this interface and add themselves; at which point users
# can configure which backend to use via settings.   Keys are the name of the library,
# which doubles as the backend= setting, and values are importable class names.
#
# NOTE: this is a temporary interface.  We DO NOT promise to continue supporting it!
#       (but if you want to experiment and don't mind breakage, here you go)
AVAILABLE_PROVIDERS = {
    "hypothesis": "hypothesis.internal.conjecture.data.HypothesisProvider",
}


class ConjectureData:
    @classmethod
    def for_buffer(
        cls,
        buffer: Union[list[int], bytes],
        *,
        observer: Optional[DataObserver] = None,
        provider: Union[type, PrimitiveProvider] = HypothesisProvider,
    ) -> "ConjectureData":
        return cls(
            len(buffer), buffer, random=None, observer=observer, provider=provider
        )

    @classmethod
    def for_choices(
        cls,
        choices: Sequence[Union[NodeTemplate, ChoiceT]],
        *,
        observer: Optional[DataObserver] = None,
        provider: Union[type, PrimitiveProvider] = HypothesisProvider,
        random: Optional[Random] = None,
    ) -> "ConjectureData":
        from hypothesis.internal.conjecture.engine import BUFFER_SIZE, choice_count

        return cls(
            max_length=BUFFER_SIZE,
            max_length_ir=choice_count(choices),
            prefix=b"",
            random=random,
            ir_prefix=choices,
            observer=observer,
            provider=provider,
        )

    def __init__(
        self,
        max_length: int,
        prefix: Union[list[int], bytes, bytearray],
        *,
        random: Optional[Random],
        observer: Optional[DataObserver] = None,
        provider: Union[type, PrimitiveProvider] = HypothesisProvider,
        ir_prefix: Optional[Sequence[Union[NodeTemplate, ChoiceT]]] = None,
        max_length_ir: Optional[int] = None,
        provider_kw: Optional[dict[str, Any]] = None,
    ) -> None:
        from hypothesis.internal.conjecture.engine import BUFFER_SIZE_IR

        if observer is None:
            observer = DataObserver()
        if provider_kw is None:
            provider_kw = {}
        elif not isinstance(provider, type):
            raise InvalidArgument(
                f"Expected {provider=} to be a class since {provider_kw=} was "
                "passed, but got an instance instead."
            )

        assert isinstance(observer, DataObserver)
        self._bytes_drawn = 0
        self.observer = observer
        self.max_length = max_length
        self.max_choices = max_length_ir
        self.max_length_ir = BUFFER_SIZE_IR
        self.is_find = False
        self.overdraw = 0
        self.__prefix = bytes(prefix)
        self.__random = random

        self.buffer: "Union[bytes, bytearray]" = bytearray()
        self.index = 0
        self.length_ir = 0
        self.index_ir = 0
        self.output = ""
        self.status = Status.VALID
        self.frozen = False
        global global_test_counter
        self.testcounter = global_test_counter
        global_test_counter += 1
        self.start_time = time.perf_counter()
        self.gc_start_time = gc_cumulative_time()
        self.events: dict[str, Union[str, int, float]] = {}
        self.interesting_origin: Optional[InterestingOrigin] = None
        self.draw_times: "dict[str, float]" = {}
        self._stateful_run_times: "defaultdict[str, float]" = defaultdict(float)
        self.max_depth = 0
        self.has_discards = False

        self.provider: PrimitiveProvider = (
            provider(self, **provider_kw) if isinstance(provider, type) else provider
        )
        assert isinstance(self.provider, PrimitiveProvider)
        if ir_prefix is None and isinstance(self.provider, HypothesisProvider):
            assert random is not None or max_length <= len(prefix)

        self.__result: "Optional[ConjectureResult]" = None

        # Observations used for targeted search.  They'll be aggregated in
        # ConjectureRunner.generate_new_examples and fed to TargetSelector.
        self.target_observations: TargetObservations = {}

        # Tags which indicate something about which part of the search space
        # this example is in. These are used to guide generation.
        self.tags: "set[StructuralCoverageTag]" = set()
        self.labels_for_structure_stack: "list[set[int]]" = []

        # Normally unpopulated but we need this in the niche case
        # that self.as_result() is Overrun but we still want the
        # examples for reporting purposes.
        self.__examples: "Optional[Examples]" = None

        # We want the top level example to have depth 0, so we start
        # at -1.
        self.depth = -1
        self.__example_record = ExampleRecord()

        # Slice indices for discrete reportable parts that which-parts-matter can
        # try varying, to report if the minimal example always fails anyway.
        self.arg_slices: set[tuple[int, int]] = set()
        self.slice_comments: dict[tuple[int, int], str] = {}
        self._observability_args: dict[str, Any] = {}
        self._observability_predicates: defaultdict = defaultdict(
            lambda: {"satisfied": 0, "unsatisfied": 0}
        )

        self.extra_information = ExtraInformation()

        self.ir_prefix = ir_prefix
        self.ir_nodes: tuple[IRNode, ...] = ()
        self.misaligned_at: Optional[MisalignedAt] = None
        self.start_example(TOP_LABEL)

    def __repr__(self) -> str:
        return "ConjectureData(%s, %d bytes%s)" % (
            self.status.name,
            len(self.buffer),
            ", frozen" if self.frozen else "",
        )

    @property
    def choices(self) -> tuple[ChoiceT, ...]:
        return tuple(node.value for node in self.ir_nodes)

    # A bit of explanation of the `observe` and `fake_forced` arguments in our
    # draw_* functions.
    #
    # There are two types of draws: sub-ir and super-ir. For instance, some ir
    # nodes use `many`, which in turn calls draw_boolean. But some strategies
    # also use many, at the super-ir level. We don't want to write sub-ir draws
    # to the DataTree (and consequently use them when computing novel prefixes),
    # since they are fully recorded by writing the ir node itself.
    # But super-ir draws are not included in the ir node, so we do want to write
    # these to the tree.
    #
    # `observe` formalizes this distinction. The draw will only be written to
    # the DataTree if observe is True.
    #
    # `fake_forced` deals with a different problem. We use `forced=` to convert
    # ir prefixes, which are potentially from other backends, into our backing
    # bits representation. This works fine, except using `forced=` in this way
    # also sets `was_forced=True` for all blocks, even those that weren't forced
    # in the traditional way. The shrinker chokes on this due to thinking that
    # nothing can be modified.
    #
    # Setting `fake_forced` to true says that yes, we want to force a particular
    # value to be returned, but we don't want to treat that block as fixed for
    # e.g. the shrinker.

    def _draw(self, ir_type, kwargs, *, observe, forced, fake_forced):
        # this is somewhat redundant with the length > max_length check at the
        # end of the function, but avoids trying to use a null self.random when
        # drawing past the node of a ConjectureData.for_choices data.
        if self.length_ir == self.max_length_ir:
            debug_report(f"overrun because hit {self.max_length_ir=}")
            self.mark_overrun()
        if len(self.ir_nodes) == self.max_choices:
            debug_report(f"overrun because hit {self.max_choices=}")
            self.mark_overrun()

        if (
            observe
            and self.ir_prefix is not None
            and self.index_ir < len(self.ir_prefix)
        ):
            choice = self._pop_choice(ir_type, kwargs, forced=forced)

            if forced is None:
                forced = choice
                fake_forced = True

        value = getattr(self.provider, f"draw_{ir_type}")(
            **kwargs, forced=forced, fake_forced=fake_forced
        )

        if observe:
            was_forced = forced is not None and not fake_forced
            getattr(self.observer, f"draw_{ir_type}")(
                value, kwargs=kwargs, was_forced=was_forced
            )
            size = 0 if self.provider.avoid_realization else ir_size([value])
            if self.length_ir + size > self.max_length_ir:
                debug_report(
                    f"overrun because {self.length_ir=} + {size=} > {self.max_length_ir=}"
                )
                self.mark_overrun()

            node = IRNode(
                ir_type=ir_type,
                value=value,
                kwargs=kwargs,
                was_forced=was_forced,
                index=len(self.ir_nodes),
            )
            self.__example_record.record_ir_draw()
            self.ir_nodes += (node,)
            self.length_ir += size

        return value

    def draw_integer(
        self,
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        weights: Optional[dict[int, float]] = None,
        shrink_towards: int = 0,
        forced: Optional[int] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> int:
        # Validate arguments
        if weights is not None:
            assert min_value is not None
            assert max_value is not None
            assert len(weights) <= 255  # arbitrary practical limit
            # We can and should eventually support total weights. But this
            # complicates shrinking as we can no longer assume we can force
            # a value to the unmapped probability mass if that mass might be 0.
            assert sum(weights.values()) < 1
            # similarly, things get simpler if we assume every value is possible.
            # we'll want to drop this restriction eventually.
            assert all(w != 0 for w in weights.values())

        if forced is not None and (min_value is None or max_value is None):
            # We draw `forced=forced - shrink_towards` here internally, after clamping.
            # If that grows larger than a 128 bit signed integer, we can't represent it.
            # Disallow this combination for now.
            # Note that bit_length() = 128 -> signed bit size = 129.
            _shrink_towards = shrink_towards
            if min_value is not None:
                _shrink_towards = max(min_value, _shrink_towards)
            if max_value is not None:
                _shrink_towards = min(max_value, _shrink_towards)

            assert (forced - _shrink_towards).bit_length() < 128
        if forced is not None and min_value is not None:
            assert min_value <= forced
        if forced is not None and max_value is not None:
            assert forced <= max_value

        kwargs: IntegerKWargs = self._pooled_kwargs(
            "integer",
            {
                "min_value": min_value,
                "max_value": max_value,
                "weights": weights,
                "shrink_towards": shrink_towards,
            },
        )
        return self._draw(
            "integer", kwargs, observe=observe, forced=forced, fake_forced=fake_forced
        )

    def draw_float(
        self,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        *,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float = SMALLEST_SUBNORMAL,
        # TODO: consider supporting these float widths at the IR level in the
        # future.
        # width: Literal[16, 32, 64] = 64,
        # exclude_min and exclude_max handled higher up,
        forced: Optional[float] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> float:
        assert smallest_nonzero_magnitude > 0
        assert not math.isnan(min_value)
        assert not math.isnan(max_value)

        if forced is not None:
            assert allow_nan or not math.isnan(forced)
            assert math.isnan(forced) or (
                sign_aware_lte(min_value, forced) and sign_aware_lte(forced, max_value)
            )

        kwargs: FloatKWargs = self._pooled_kwargs(
            "float",
            {
                "min_value": min_value,
                "max_value": max_value,
                "allow_nan": allow_nan,
                "smallest_nonzero_magnitude": smallest_nonzero_magnitude,
            },
        )
        return self._draw(
            "float", kwargs, observe=observe, forced=forced, fake_forced=fake_forced
        )

    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        forced: Optional[str] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> str:
        assert forced is None or min_size <= len(forced) <= max_size
        assert min_size >= 0

        kwargs: StringKWargs = self._pooled_kwargs(
            "string",
            {
                "intervals": intervals,
                "min_size": min_size,
                "max_size": max_size,
            },
        )
        return self._draw(
            "string", kwargs, observe=observe, forced=forced, fake_forced=fake_forced
        )

    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        *,
        forced: Optional[bytes] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> bytes:
        assert forced is None or min_size <= len(forced) <= max_size
        assert min_size >= 0

        kwargs: BytesKWargs = self._pooled_kwargs(
            "bytes", {"min_size": min_size, "max_size": max_size}
        )
        return self._draw(
            "bytes", kwargs, observe=observe, forced=forced, fake_forced=fake_forced
        )

    def draw_boolean(
        self,
        p: float = 0.5,
        *,
        forced: Optional[bool] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> bool:
        # Internally, we treat probabilities lower than 1 / 2**64 as
        # unconditionally false.
        #
        # Note that even if we lift this 64 bit restriction in the future, p
        # cannot be 0 (1) when forced is True (False).
        eps = 2 ** (-64) if isinstance(self.provider, HypothesisProvider) else 0
        assert (forced is not True) or (0 + eps) < p
        assert (forced is not False) or p < (1 - eps)

        kwargs: BooleanKWargs = self._pooled_kwargs("boolean", {"p": p})
        return self._draw(
            "boolean", kwargs, observe=observe, forced=forced, fake_forced=fake_forced
        )

    def _pooled_kwargs(self, ir_type, kwargs):
        """Memoize common dictionary objects to reduce memory pressure."""
        # caching runs afoul of nondeterminism checks
        if self.provider.avoid_realization:
            return kwargs

        key = (ir_type, *choice_kwargs_key(ir_type, kwargs))
        try:
            return POOLED_KWARGS_CACHE[key]
        except KeyError:
            POOLED_KWARGS_CACHE[key] = kwargs
            return kwargs

    def _pop_choice(
        self, ir_type: ChoiceNameT, kwargs: ChoiceKwargsT, *, forced: Optional[ChoiceT]
    ) -> ChoiceT:
        assert self.ir_prefix is not None
        # checked in _draw
        assert self.index_ir < len(self.ir_prefix)

        value = self.ir_prefix[self.index_ir]
        if isinstance(value, NodeTemplate):
            node: NodeTemplate = value
            if node.count is not None:
                assert node.count >= 0
            # node templates have to be at the end for now, since it's not immediately
            # apparent how to handle overruning a node template while generating a single
            # node if the alternative is not "the entire data is an overrun".
            assert self.index_ir == len(self.ir_prefix) - 1
            if node.type == "simplest":
                if isinstance(self.provider, HypothesisProvider):
                    try:
                        choice: ChoiceT = choice_from_index(0, ir_type, kwargs)
                    except ChoiceTooLarge:
                        self.mark_overrun()
                else:
                    # give alternative backends control over these draws
                    choice = getattr(self.provider, f"draw_{ir_type}")(**kwargs)
            else:
                raise NotImplementedError

            if node.count is not None:
                node.count -= 1
                if node.count < 0:
                    self.mark_overrun()
            return choice

        choice = value
        node_ir_type = {
            str: "string",
            float: "float",
            int: "integer",
            bool: "boolean",
            bytes: "bytes",
        }[type(choice)]
        # If we're trying to:
        # * draw a different ir type at the same location
        # * draw the same ir type with a different kwargs, which does not permit
        #   the current value
        #
        # then we call this a misalignment, because the choice sequence has
        # slipped from what we expected at some point. An easy misalignment is
        #
        #   one_of(integers(0, 100), integers(101, 200))
        #
        # where the choice sequence [0, 100] has kwargs {min_value: 0, max_value: 100}
        # at index 1, but [0, 101] has kwargs {min_value: 101, max_value: 200} at
        # index 1 (which does not permit any of the values 0-100).
        #
        # When the choice sequence becomes misaligned, we generate a new value of the
        # type and kwargs the strategy expects.
        if node_ir_type != ir_type or not choice_permitted(choice, kwargs):
            # only track first misalignment for now.
            if self.misaligned_at is None:
                self.misaligned_at = (self.index_ir, ir_type, kwargs, forced)
            try:
                # Fill in any misalignments with index 0 choices. An alternative to
                # this is using the index of the misaligned choice instead
                # of index 0, which may be useful for maintaining
                # "similarly-complex choices" in the shrinker. This requires
                # attaching an index to every choice in ConjectureData.for_choices,
                # which we don't always have (e.g. when reading from db).
                #
                # If we really wanted this in the future we could make this complexity
                # optional, use it if present, and default to index 0 otherwise.
                # This complicates our internal api and so I'd like to avoid it
                # if possible.
                #
                # Additionally, I don't think slips which require
                # slipping to high-complexity values are common. Though arguably
                # we may want to expand a bit beyond *just* the simplest choice.
                # (we could for example consider sampling choices from index 0-10).
                choice = choice_from_index(0, ir_type, kwargs)
            except ChoiceTooLarge:
                # should really never happen with a 0-index choice, but let's be safe.
                self.mark_overrun()

        self.index_ir += 1
        return choice

    def as_result(self) -> Union[ConjectureResult, _Overrun]:
        """Convert the result of running this test into
        either an Overrun object or a ConjectureResult."""

        assert self.frozen
        if self.status == Status.OVERRUN:
            return Overrun
        if self.__result is None:
            self.__result = ConjectureResult(
                status=self.status,
                interesting_origin=self.interesting_origin,
                buffer=self.buffer,
                examples=self.examples,
                ir_nodes=self.ir_nodes,
                output=self.output,
                extra_information=(
                    self.extra_information
                    if self.extra_information.has_information()
                    else None
                ),
                has_discards=self.has_discards,
                target_observations=self.target_observations,
                tags=frozenset(self.tags),
                arg_slices=self.arg_slices,
                slice_comments=self.slice_comments,
                misaligned_at=self.misaligned_at,
            )
            assert self.__result is not None
        return self.__result

    def __assert_not_frozen(self, name: str) -> None:
        if self.frozen:
            raise Frozen(f"Cannot call {name} on frozen ConjectureData")

    def note(self, value: Any) -> None:
        self.__assert_not_frozen("note")
        if not isinstance(value, str):
            value = repr(value)
        self.output += value

    def draw(
        self,
        strategy: "SearchStrategy[Ex]",
        label: Optional[int] = None,
        observe_as: Optional[str] = None,
    ) -> "Ex":
        from hypothesis.internal.observability import TESTCASE_CALLBACKS
        from hypothesis.strategies._internal.utils import to_jsonable

        if self.is_find and not strategy.supports_find:
            raise InvalidArgument(
                f"Cannot use strategy {strategy!r} within a call to find "
                "(presumably because it would be invalid after the call had ended)."
            )

        at_top_level = self.depth == 0
        start_time = None
        if at_top_level:
            # We start this timer early, because accessing attributes on a LazyStrategy
            # can be almost arbitrarily slow.  In cases like characters() and text()
            # where we cache something expensive, this led to Flaky deadline errors!
            # See https://github.com/HypothesisWorks/hypothesis/issues/2108
            start_time = time.perf_counter()
            gc_start_time = gc_cumulative_time()

        strategy.validate()

        if strategy.is_empty:
            self.mark_invalid(f"empty strategy {self!r}")

        if self.depth >= MAX_DEPTH:
            self.mark_invalid("max depth exceeded")

        if label is None:
            assert isinstance(strategy.label, int)
            label = strategy.label
        self.start_example(label=label)
        try:
            if not at_top_level:
                return strategy.do_draw(self)
            assert start_time is not None
            key = observe_as or f"generate:unlabeled_{len(self.draw_times)}"
            try:
                strategy.validate()
                try:
                    v = strategy.do_draw(self)
                finally:
                    # Subtract the time spent in GC to avoid overcounting, as it is
                    # accounted for at the overall example level.
                    in_gctime = gc_cumulative_time() - gc_start_time
                    self.draw_times[key] = time.perf_counter() - start_time - in_gctime
            except Exception as err:
                add_note(
                    err,
                    f"while generating {key.removeprefix('generate:')!r} from {strategy!r}",
                )
                raise
            if TESTCASE_CALLBACKS:
                self._observability_args[key] = to_jsonable(v)
            return v
        finally:
            self.stop_example()

    def start_example(self, label: int) -> None:
        self.provider.span_start(label)
        self.__assert_not_frozen("start_example")
        self.depth += 1
        # Logically it would make sense for this to just be
        # ``self.depth = max(self.depth, self.max_depth)``, which is what it used to
        # be until we ran the code under tracemalloc and found a rather significant
        # chunk of allocation was happening here. This was presumably due to varargs
        # or the like, but we didn't investigate further given that it was easy
        # to fix with this check.
        if self.depth > self.max_depth:
            self.max_depth = self.depth
        self.__example_record.start_example(label)
        self.labels_for_structure_stack.append({label})

    def stop_example(self, *, discard: bool = False) -> None:
        self.provider.span_end(discard)
        if self.frozen:
            return
        if discard:
            self.has_discards = True
        self.depth -= 1
        assert self.depth >= -1
        self.__example_record.stop_example(discard=discard)

        labels_for_structure = self.labels_for_structure_stack.pop()

        if not discard:
            if self.labels_for_structure_stack:
                self.labels_for_structure_stack[-1].update(labels_for_structure)
            else:
                self.tags.update([structural_coverage(l) for l in labels_for_structure])

        if discard:
            # Once we've discarded an example, every test case starting with
            # this prefix contains discards. We prune the tree at that point so
            # as to avoid future test cases bothering with this region, on the
            # assumption that some example that you could have used instead
            # there would *not* trigger the discard. This greatly speeds up
            # test case generation in some cases, because it allows us to
            # ignore large swathes of the search space that are effectively
            # redundant.
            #
            # A scenario that can cause us problems but which we deliberately
            # have decided not to support is that if there are side effects
            # during data generation then you may end up with a scenario where
            # every good test case generates a discard because the discarded
            # section sets up important things for later. This is not terribly
            # likely and all that you see in this case is some degradation in
            # quality of testing, so we don't worry about it.
            #
            # Note that killing the branch does *not* mean we will never
            # explore below this point, and in particular we may do so during
            # shrinking. Any explicit request for a data object that starts
            # with the branch here will work just fine, but novel prefix
            # generation will avoid it, and we can use it to detect when we
            # have explored the entire tree (up to redundancy).

            self.observer.kill_branch()

    @property
    def examples(self) -> Examples:
        assert self.frozen
        if self.__examples is None:
            self.__examples = Examples(record=self.__example_record)
        return self.__examples

    def freeze(self) -> None:
        if self.frozen:
            assert isinstance(self.buffer, bytes)
            return
        self.finish_time = time.perf_counter()
        self.gc_finish_time = gc_cumulative_time()
        assert len(self.buffer) == self.index

        # Always finish by closing all remaining examples so that we have a
        # valid tree.
        while self.depth >= 0:
            self.stop_example()

        self.__example_record.freeze()
        self.frozen = True
        self.buffer = bytes(self.buffer)
        self.observer.conclude_test(self.status, self.interesting_origin)

    def choice(
        self,
        values: Sequence[T],
        *,
        forced: Optional[T] = None,
        fake_forced: bool = False,
        observe: bool = True,
    ) -> T:
        forced_i = None if forced is None else values.index(forced)
        i = self.draw_integer(
            0,
            len(values) - 1,
            forced=forced_i,
            fake_forced=fake_forced,
            observe=observe,
        )
        return values[i]

    def draw_bits(
        self, n: int, *, forced: Optional[int] = None, fake_forced: bool = False
    ) -> int:
        """Return an ``n``-bit integer from the underlying source of
        bytes. If ``forced`` is set to an integer will instead
        ignore the underlying source and simulate a draw as if it had
        returned that integer."""
        self.__assert_not_frozen("draw_bits")
        if n == 0:
            return 0
        assert n > 0
        n_bytes = bits_to_bytes(n)
        self.__check_capacity(n_bytes)

        if forced is not None:
            buf = int_to_bytes(forced, n_bytes)
        elif self._bytes_drawn < len(self.__prefix):
            index = self._bytes_drawn
            buf = self.__prefix[index : index + n_bytes]
            if len(buf) < n_bytes:  # pragma: no cover # removing soon
                assert self.__random is not None
                buf += uniform(self.__random, n_bytes - len(buf))
        else:
            assert self.__random is not None
            buf = uniform(self.__random, n_bytes)
        buf = bytearray(buf)
        self._bytes_drawn += n_bytes

        assert len(buf) == n_bytes

        # If we have a number of bits that is not a multiple of 8
        # we have to mask off the high bits.
        buf[0] &= BYTE_MASKS[n % 8]
        buf = bytes(buf)
        result = int_from_bytes(buf)

        assert isinstance(self.buffer, bytearray)
        self.buffer.extend(buf)
        self.index = len(self.buffer)

        assert result.bit_length() <= n
        return result

    def __check_capacity(self, n: int) -> None:
        if self.index + n > self.max_length:
            self.mark_overrun()

    def conclude_test(
        self,
        status: Status,
        interesting_origin: Optional[InterestingOrigin] = None,
    ) -> NoReturn:
        assert (interesting_origin is None) or (status == Status.INTERESTING)
        self.__assert_not_frozen("conclude_test")
        self.interesting_origin = interesting_origin
        self.status = status
        self.freeze()
        raise StopTest(self.testcounter)

    def mark_interesting(
        self, interesting_origin: Optional[InterestingOrigin] = None
    ) -> NoReturn:
        self.conclude_test(Status.INTERESTING, interesting_origin)

    def mark_invalid(self, why: Optional[str] = None) -> NoReturn:
        if why is not None:
            self.events["invalid because"] = why
        self.conclude_test(Status.INVALID)

    def mark_overrun(self) -> NoReturn:
        self.conclude_test(Status.OVERRUN)


def bits_to_bytes(n: int) -> int:
    """The number of bytes required to represent an n-bit number.
    Equivalent to (n + 7) // 8, but slightly faster. This really is
    called enough times that that matters."""
    return (n + 7) >> 3


def draw_choice(ir_type, kwargs, *, random):
    from hypothesis.internal.conjecture.engine import BUFFER_SIZE

    cd = ConjectureData(
        max_length=BUFFER_SIZE,
        # buffer doesn't matter if forced is passed since we're forcing the sole draw
        prefix=b"",
        random=random,
    )
    return getattr(cd.provider, f"draw_{ir_type}")(**kwargs)
