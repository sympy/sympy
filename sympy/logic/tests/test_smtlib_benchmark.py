"""Tests for the SMT-LIB benchmarking utility.

The utility is meant to be pointed at a benchmark suite, which is far too
large to ship or to run in CI; these tests only check that it reports the right
thing for each kind of input.
"""

from __future__ import annotations

import os
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from sympy.external import import_module
from sympy.logic.utilities.smtlib_benchmark import (
    FileResult, collect_files, main, run_file, run_file_with_timeout
)

lark = import_module('lark')

# disable tests if lark is not present
disabled = lark is None


SAT_LRA = """
(set-logic QF_LRA)
(set-info :status sat)
(declare-const x Real)
(declare-const y Real)
(assert (> (+ x y) 3))
(assert (< x 1))
(check-sat)
(exit)
"""

UNSAT_LRA = """
(set-logic QF_LRA)
(set-info :status unsat)
(declare-const x Real)
(assert (> x 3))
(assert (< x 1))
(check-sat)
"""

UNSAT_BOOL = """
(set-logic QF_UF)
(set-info :status unsat)
(declare-const p Bool)
(assert (and p (not p)))
(check-sat)
"""

MISLABELLED = """
(set-logic QF_LRA)
(set-info :status unsat)
(declare-const x Real)
(assert (> x 3))
(check-sat)
"""

NO_STATUS = """
(set-logic QF_LRA)
(declare-const x Real)
(assert (> x 3))
(check-sat)
"""

QUANTIFIED = """
(set-logic LRA)
(declare-fun x () Real)
(assert (forall ((y Real)) (> y x)))
(check-sat)
"""

BITVECTOR = """
(set-logic QF_BV)
(declare-const a (_ BitVec 8))
(assert (= (bvadd a a) a))
(check-sat)
"""

MIXED_BOOL_AND_THEORY = """
(set-logic QF_LRA)
(declare-const p Bool)
(declare-const x Real)
(assert (or p (> x 3)))
(assert (< x 1))
(check-sat)
"""

MALFORMED = """
(declare-const x Real
(assert (> x 3))
"""


def _write(directory, name, source):
    path = Path(directory) / name
    path.write_text(source, encoding='utf-8')
    return str(path)


def test_sat_and_unsat_are_checked_against_the_recorded_status():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'sat.smt2', SAT_LRA))
        assert result.status == 'sat'
        assert result.expected == 'sat'
        assert result.check == 'ok'

        result = run_file(_write(directory, 'unsat.smt2', UNSAT_LRA))
        assert result.status == 'unsat'
        assert result.check == 'ok'


def test_a_wrong_answer_is_reported():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'bad.smt2', MISLABELLED))
        assert result.status == 'sat'
        assert result.expected == 'unsat'
        assert result.check == 'WRONG'


def test_a_file_without_a_status_is_not_checked():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'plain.smt2', NO_STATUS))
        assert result.status == 'sat'
        assert result.expected == ''
        assert result.check == ''


def test_purely_boolean_problems_are_solved():
    # the LRA theory solver cannot take propositional variables, so a problem
    # without theory atoms has to fall back to plain dpll2
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'bool.smt2', UNSAT_BOOL))
        assert result.status == 'unsat'
        assert result.check == 'ok'


def test_every_phase_is_timed():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'sat.smt2', SAT_LRA))
        assert set(result.times) == {'parse', 'build', 'solve'}
        assert all(t >= 0 for t in result.times.values())
        assert result.total >= sum(result.times.values())


def test_phases_that_did_not_run_are_not_timed():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'broken.smt2', MALFORMED))
        assert result.status == 'parse-error'
        assert 'solve' not in result.times


def test_solver_none_only_parses():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'sat.smt2', SAT_LRA), solver='none')
        assert result.status == 'parsed'
        assert result.check == ''
        assert 'solve' not in result.times


def test_unsupported_input_is_reported_with_a_reason():
    with tempfile.TemporaryDirectory() as directory:
        for name, source, expected in [('quant.smt2', QUANTIFIED, 'Quantifiers'),
                                       ('bv.smt2', BITVECTOR, 'bvadd')]:
            result = run_file(_write(directory, name, source))
            assert result.status == 'unsupported', name
            assert expected in result.reason, name


def test_malformed_input_is_reported_as_a_parse_error():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'broken.smt2', MALFORMED))
        assert result.status == 'parse-error'
        assert 'line 3' in result.reason
        # the reason has to stay on one line so that the report lines up
        assert '\n' not in result.reason


def test_an_unexpected_failure_is_reported_rather_than_raised():
    # the LRA solver raises ValueError on a problem that mixes Boolean
    # variables with theory atoms; a whole run must not stop because of it
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(_write(directory, 'mixed.smt2', MIXED_BOOL_AND_THEORY))
        assert result.status == 'error'
        assert 'Unhandled Predicate' in result.reason


def test_missing_file_is_reported_as_an_error():
    with tempfile.TemporaryDirectory() as directory:
        result = run_file(os.path.join(directory, 'nope.smt2'))
        assert result.status == 'error'
        assert 'could not read file' in result.reason


def test_collect_files_searches_directories_recursively():
    with tempfile.TemporaryDirectory() as directory:
        _write(directory, 'a.smt2', SAT_LRA)
        _write(directory, 'notes.txt', 'not a benchmark')
        nested = Path(directory) / 'nested'
        nested.mkdir()
        _write(nested, 'b.smt', UNSAT_LRA)

        found = [os.path.basename(p) for p in collect_files([directory])]
        assert sorted(found) == ['a.smt2', 'b.smt']

        # a file given by name is used whatever its extension
        assert len(collect_files([os.path.join(directory, 'notes.txt')])) == 1


def test_a_worker_that_runs_out_of_time_reports_a_timeout():
    with tempfile.TemporaryDirectory() as directory:
        path = _write(directory, 'sat.smt2', SAT_LRA)
        result = run_file_with_timeout(path, timeout=0)
        assert result.status == 'timeout'
        assert 'gave up' in result.reason


def test_a_file_solved_in_a_worker_gives_the_same_answer():
    with tempfile.TemporaryDirectory() as directory:
        path = _write(directory, 'unsat.smt2', UNSAT_LRA)
        result = run_file_with_timeout(path, timeout=60)
        assert result.status == 'unsat'
        assert result.check == 'ok'
        assert result.times['solve'] > 0


def test_check_needs_an_answer_and_an_expected_status():
    assert FileResult(path='a', status='sat', expected='sat').check == 'ok'
    assert FileResult(path='a', status='sat', expected='unsat').check == 'WRONG'
    assert FileResult(path='a', status='sat').check == ''
    assert FileResult(path='a', status='timeout', expected='sat').check == ''
    assert FileResult(path='a', status='unknown', expected='sat').check == ''


def _main(argv):
    """Run the command line entry point without printing its report."""
    with open(os.devnull, 'w') as sink:
        with redirect_stdout(sink), redirect_stderr(sink):
            return main(argv)


def test_main_returns_nonzero_only_for_a_wrong_answer():
    with tempfile.TemporaryDirectory() as directory:
        assert _main([_write(directory, 'sat.smt2', SAT_LRA)]) == 0
        assert _main([_write(directory, 'bad.smt2', MISLABELLED)]) == 1
        # files that cannot be handled are not failures
        assert _main([_write(directory, 'q.smt2', QUANTIFIED)]) == 0


def test_main_rejects_a_directory_with_no_benchmarks():
    with tempfile.TemporaryDirectory() as directory:
        assert _main([directory]) == 2
