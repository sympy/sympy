"""Command line utility for benchmarking SymPy's SAT/SMT solvers on SMT-LIB
files.

There are standard benchmark suites for SMT solvers (see https://smt-lib.org).
This runs SymPy's solvers over such a suite, checks the answers against the
``(set-info :status ...)`` line recorded in each file and reports how long
parsing, building the SymPy expression and solving each took::

    $ python -m sympy.logic.utilities.smtlib_benchmark --solver lra -t 5 QF_LRA/

Each ``PATH`` is either an SMT-LIB file or a directory that is searched
recursively for them. Files are run in a child process so that a solver taking
too long can be killed; call :func:`run_file` directly to run one in process,
which is easier to debug and to profile.

This is a developer tool rather than part of the test suite: benchmark files
are not distributed with SymPy and running a suite takes far longer than a test
run.

Note that SymPy has no solver for integer arithmetic, so a problem declaring
``Int`` constants is solved over the reals and may well disagree with the
answer the file records.
"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import queue
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from sympy.assumptions.ask import Q
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.logic.algorithms.lra_theory import UnhandledInput
from sympy.logic.algorithms.z3_wrapper import z3_satisfiable
from sympy.logic.boolalg import And, BooleanAtom, BooleanFunction
from sympy.logic.inference import satisfiable
from sympy.parsing.smtlib.lark.smtlib_parser import (
    LarkSMTLibParser, SMTLibSyntaxError, UnknownSMTLibCommandError,
    UnknownSMTLibOperatorError)
from sympy.parsing.smtlib.lark.transformer import SMTLibTransformer


#: Solvers that ``--solver`` accepts. ``none`` only parses the files, which is
#: useful for measuring the parser on its own.
SOLVERS = ('lra', 'z3', 'dpll2', 'dpll', 'pycosat', 'minisat22', 'none')

#: Filename extensions searched for when a directory is given.
EXTENSIONS = ('.smt2', '.smt')

#: The steps that are timed separately, in the order they run.
PHASES = ('parse', 'build', 'solve')

#: Modules the solvers that are not built into SymPy need.
SOLVER_MODULES = {'z3': 'z3', 'pycosat': 'pycosat', 'minisat22': 'pysat'}

DEFAULT_TIMEOUT = 30.0


@dataclass
class FileResult:
    """The outcome of running one SMT-LIB file.

    Every field is a plain Python object so that a result can be sent back
    from the child process that produced it.
    """

    path: str
    #: ``sat``, ``unsat``, ``unknown``, ``parsed``, ``timeout``,
    #: ``unsupported``, ``parse-error`` or ``error``
    status: str = 'error'
    #: the answer recorded by ``(set-info :status ...)``, or ``''``
    expected: str = ''
    #: why the file could not be handled, for the non-answer statuses
    reason: str = ''
    times: dict[str, float] = field(default_factory=dict)
    total: float = 0.0

    @property
    def check(self):
        """``ok`` or ``WRONG`` when the file records an expected answer."""
        if self.status not in ('sat', 'unsat') or self.expected not in ('sat', 'unsat'):
            return ''
        return 'ok' if self.status == self.expected else 'WRONG'


def collect_files(paths):
    """Expand ``paths`` into a sorted list of SMT-LIB files.

    Directories are searched recursively; a file named directly is used
    whatever its extension.
    """
    files = []
    for path in map(Path, paths):
        if path.is_dir():
            files.extend(sorted(p for p in path.rglob('*')
                                if p.is_file() and p.suffix.lower() in EXTENSIONS))
        else:
            files.append(path)
    return files


_parser = None


def get_parser():
    """Return the shared SMT-LIB parser, building it on first use.

    Building it compiles the grammar, which is slow enough to be worth doing
    once up front rather than charging it to the first file.
    """
    global _parser
    if _parser is None:
        _parser = LarkSMTLibParser(transform=False)
    return _parser


def run_file(path, solver='lra'):
    """Parse and solve one SMT-LIB file and return a :class:`FileResult`.

    This runs in the calling process and has no timeout of its own; use
    :func:`run_file_with_timeout` to get one.

    ``solver`` is one of :data:`SOLVERS`. ``none`` parses the file without
    solving it, and ``lra`` means dpll2 with the linear arithmetic theory,
    falling back to plain dpll2 for problems that have no theory atoms.

    Whatever a stage raises becomes the file's reported outcome rather than
    propagating: finding out how SymPy fails on a file is the point of the
    exercise, and one unhandled file must not stop a run over a whole suite.
    That is why the handlers below catch ``Exception`` broadly.
    """
    result = FileResult(path=str(path))
    started = time.perf_counter()

    def finish(status, reason=''):
        result.status = status
        result.reason = reason
        result.total = time.perf_counter() - started
        return result

    def timed(name, func):
        start = time.perf_counter()
        try:
            return func()
        finally:
            result.times[name] = time.perf_counter() - start

    try:
        text = Path(path).read_text(encoding='utf-8', errors='replace')
    except OSError as exc:
        return finish('error', f"could not read file: {exc}")

    try:
        tree = timed('parse', lambda: get_parser().doparse(text))
    except SMTLibSyntaxError as exc:
        return finish('parse-error', _describe(exc))
    except Exception as exc:  # noqa: BLE001
        return finish('error', _describe(exc))

    def build():
        transformer = SMTLibTransformer()
        transformer.transform(tree)
        return transformer

    try:
        transformer = timed('build', build)
    except Exception as raised:  # noqa: BLE001
        # lark wraps whatever the transformer raised in a VisitError
        cause = getattr(raised, 'orig_exc', raised)
        unsupported = (NotImplementedError, UnknownSMTLibCommandError,
                       UnknownSMTLibOperatorError)
        return finish('unsupported' if isinstance(cause, unsupported) else 'error',
                      _describe(cause))

    result.expected = str(transformer.get_info().get(':status') or '')
    # (declare-const x Real) is parsed as a Q.real(x) assertion, which is an
    # assumption about the symbol rather than a constraint the solvers take
    expr = And(*[a for a in transformer.get_result()[1]
                 if not (isinstance(a, AppliedPredicate)
                         and a.function in (Q.real, Q.integer))])

    if solver == 'none':
        return finish('parsed')
    if solver == 'lra' and not _has_theory_atoms(expr):
        solver = 'dpll2'  # the LRA theory solver rejects bare propositions

    try:
        model = timed('solve', lambda: _solve(expr, solver))
    except UnhandledInput as exc:
        return finish('unsupported', _describe(exc))
    except Exception as exc:  # noqa: BLE001
        return finish('error', _describe(exc))

    if model is None:
        return finish('unknown')  # z3 returns None when it gives up
    # a satisfying model can be an empty dict, so compare against False
    return finish('unsat' if model is False else 'sat')


def _has_theory_atoms(expr):
    """Whether the formula says anything beyond relating its propositions.

    An atom is either a bare symbol used as a proposition or something a theory
    has to reason about, such as ``x < 1``; only the latter needs LRA.
    """
    stack = [expr]
    while stack:
        node = stack.pop()
        if isinstance(node, BooleanFunction):
            stack.extend(node.args)
        elif not isinstance(node, (Symbol, BooleanAtom)):
            return True
    return False


def _solve(expr, solver):
    """Run ``solver`` on ``expr`` and return its raw result."""
    if solver == 'z3':
        return z3_satisfiable(expr)
    if solver == 'lra':
        return satisfiable(expr, use_lra_theory=True)
    return satisfiable(expr, algorithm=solver)


#: How much of a reason to show. z3 reports every sort error in the problem at
#: once, which runs to thousands of characters.
MAX_REASON = 200


def _describe(exc):
    """A one-line description of ``exc``, so the report stays lined up."""
    message = ' '.join(str(exc).split())
    if message:
        text = f"{type(exc).__name__}: {message}"
    else:
        # a bare assert has nothing to say for itself, so say where it came from
        frame = exc.__traceback__
        while frame is not None and frame.tb_next is not None:
            frame = frame.tb_next
        if frame is None:
            return type(exc).__name__
        code = frame.tb_frame.f_code
        text = (f"{type(exc).__name__} in {code.co_name} at "
                f"{os.path.basename(code.co_filename)}:{frame.tb_lineno}")
    return text if len(text) <= MAX_REASON else text[:MAX_REASON - 3] + '...'


def _worker(path, solver, results):
    results.put(run_file(path, solver))


def run_file_with_timeout(path, solver='lra', timeout=DEFAULT_TIMEOUT):
    """Run one file in a child process, killing it after ``timeout`` seconds.

    Forking, where it is available, lets the child inherit the compiled grammar
    and the already imported SymPy, which otherwise cost more than solving most
    files.
    """
    try:
        ctx = multiprocessing.get_context('fork')
    except ValueError:
        ctx = multiprocessing.get_context()

    results = ctx.Queue()
    process = ctx.Process(target=_worker, args=(str(path), solver, results))
    started = time.perf_counter()
    process.start()
    try:
        result = results.get(timeout=timeout)
    except queue.Empty:
        result = FileResult(path=str(path), status='timeout',
                            reason=f"gave up after {timeout:g}s")
    finally:
        if process.is_alive():
            process.terminate()
        process.join(5)

    result.total = time.perf_counter() - started
    return result


_ROW = '%-12s %-9s %-6s %s  %s'


def _print_row(cells, path):
    print(_ROW % (cells[0], cells[1], cells[2], ' '.join(cells[3:]), path))
    sys.stdout.flush()


def _print_result(result):
    times = ['%8.3f' % result.times[p] if p in result.times else ' ' * 8
             for p in PHASES]
    _print_row([result.status, result.expected or '-', result.check or '-',
                *times, '%8.3f' % result.total], result.path)
    if result.reason:
        print(f"    {result.reason}")


def _print_summary(results):
    # a file that was answered and checked counts as ok or WRONG; one with no
    # recorded status counts under the answer it got
    counts = Counter(r.check or r.status for r in results)
    print('\nsummary\n-------')
    print('files             %d' % len(results))
    for name, count in sorted(counts.items()):
        print('%-17s %d' % (name, count))
    phases = ', '.join('%s %.3fs' % (p, sum(r.times.get(p, 0.0) for r in results))
                       for p in PHASES)
    print('time              %.3fs (%s)' % (sum(r.total for r in results), phases))


def main(argv=None):
    """Entry point for ``python -m sympy.logic.utilities.smtlib_benchmark``."""
    parser = argparse.ArgumentParser(
        prog='python -m sympy.logic.utilities.smtlib_benchmark',
        description=__doc__.split('\n\n')[0].strip())
    parser.add_argument('paths', metavar='PATH', nargs='+',
                        help='SMT-LIB files, or directories to search recursively')
    parser.add_argument('-s', '--solver', default='lra', choices=SOLVERS,
                        help='solver to run (default: %(default)s); "none" only '
                             'parses the files')
    parser.add_argument('-t', '--timeout', type=float, default=DEFAULT_TIMEOUT,
                        metavar='SECONDS',
                        help='give up on a file after this long (default: %(default)s)')
    args = parser.parse_args(argv)

    # satisfiable() silently falls back to dpll2 when the module a solver needs
    # is missing, which would quietly benchmark the wrong solver
    for module in ('lark', SOLVER_MODULES.get(args.solver)):
        if module and import_module(module) is None:
            print(f"this needs the {module!r} module, which is not installed",
                  file=sys.stderr)
            return 2

    files = collect_files(args.paths)
    if not files:
        print('no SMT-LIB files found (looked for %s)' % ' '.join(EXTENSIONS),
              file=sys.stderr)
        return 2

    # build the parser here: every forked child inherits it, so its cost
    # belongs to the run as a whole rather than to the first file
    get_parser()
    _print_row(['status', 'expected', 'check',
                *['%8s' % p for p in PHASES], '%8s' % 'total'], 'file')

    results = []
    for path in files:
        result = run_file_with_timeout(path, args.solver, args.timeout)
        results.append(result)
        _print_result(result)

    _print_summary(results)
    return 1 if any(r.check == 'WRONG' for r in results) else 0


if __name__ == '__main__':
    sys.exit(main())
