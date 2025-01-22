# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import math
from typing import Optional

from hypothesis.internal.compat import int_from_bytes
from hypothesis.internal.conjecture.data import (
    BYTE_MASKS,
    COLLECTION_DEFAULT_MAX_SIZE,
    ConjectureData,
    PrimitiveProvider,
    bits_to_bytes,
)
from hypothesis.internal.conjecture.floats import lex_to_float
from hypothesis.internal.conjecture.utils import many
from hypothesis.internal.floats import make_float_clamper
from hypothesis.internal.intervalsets import IntervalSet


class BytestringProvider(PrimitiveProvider):
    lifetime = "test_case"

    def __init__(
        self, conjecturedata: Optional["ConjectureData"], /, *, bytestring: bytes
    ):
        super().__init__(conjecturedata)
        self.bytestring = bytestring
        self.index = 0
        self.drawn = bytearray()

    def _draw_bits(self, n):
        if n == 0:  # pragma: no cover
            return 0
        n_bytes = bits_to_bytes(n)
        if self.index + n_bytes > len(self.bytestring):
            self._cd.mark_overrun()
        buf = bytearray(self.bytestring[self.index : self.index + n_bytes])
        self.index += n_bytes

        buf[0] &= BYTE_MASKS[n % 8]
        buf = bytes(buf)
        self.drawn += buf
        return int_from_bytes(buf)

    def draw_boolean(
        self,
        p: float = 0.5,
        *,
        forced: Optional[bool] = None,
        fake_forced: bool = False,
    ) -> bool:
        if forced is not None:
            return forced

        if p <= 0:
            return False
        if p >= 1:
            return True

        # always use one byte for booleans to maintain constant draw size.
        # If a probability requires more than 8 bits to represent precisely,
        # the result will be slightly biased, but not badly.
        bits = 8
        size = 2**bits
        # always leave at least one value that can be true, even for very small
        # p.
        falsey = max(1, math.floor(size * (1 - p)))
        n = self._draw_bits(bits)
        return n >= falsey

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
        if forced is not None:
            return forced

        assert self._cd is not None

        # we explicitly ignore integer weights for now, as they are likely net
        # negative on fuzzer performance.

        if min_value is None and max_value is None:
            min_value = -(2**127)
            max_value = 2**127 - 1
        elif min_value is None:
            assert max_value is not None
            min_value = max_value - 2**64
        elif max_value is None:
            assert min_value is not None
            max_value = min_value + 2**64

        if min_value == max_value:
            return min_value

        bits = (max_value - min_value).bit_length()
        value = self._draw_bits(bits)
        while not (min_value <= value <= max_value):
            value = self._draw_bits(bits)
        return value

    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
        forced: Optional[float] = None,
        fake_forced: bool = False,
    ) -> float:
        if forced is not None:
            return forced

        n = self._draw_bits(64)
        sign = -1 if n >> 64 else 1
        f = sign * lex_to_float(n & ((1 << 64) - 1))
        clamper = make_float_clamper(
            min_value,
            max_value,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
            allow_nan=allow_nan,
        )
        return clamper(f)

    def _draw_collection(self, min_size, max_size, *, alphabet_size):
        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            observe=False,
        )
        values = []
        while elements.more():
            values.append(self.draw_integer(0, alphabet_size - 1))
        return values

    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        forced: Optional[str] = None,
        fake_forced: bool = False,
    ) -> str:
        if forced is not None:
            return forced
        values = self._draw_collection(min_size, max_size, alphabet_size=len(intervals))
        return "".join(chr(intervals[v]) for v in values)

    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
        *,
        forced: Optional[bytes] = None,
        fake_forced: bool = False,
    ) -> bytes:
        if forced is not None:
            return forced
        values = self._draw_collection(min_size, max_size, alphabet_size=2**8)
        return bytes(values)
