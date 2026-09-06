#!/usr/bin/env python
"""Run the Rubi integration test suite through manualintegrate.

For every ``RubiTestSuiteCase`` in ``rubi_rules/rubi_test_suite`` (in the
matchpy checkout next to this repository) this script records:

* the rule tree returned by ``integral_steps(integrand, variable)``,
* whether the tree contains ``DontKnowRule``,
* the antiderivative computed by ``manualintegrate(integrand, variable)``,
* whether that antiderivative is correct, verified by differentiating it
  and checking that the difference with the integrand simplifies (or
  numerically evaluates) to zero.

Results are written as JSON Lines, one file per test module, mirroring the
directory structure of the test suite under ``--results-dir``.  SymPy
objects are serialized with ``str(expr)`` whenever
``eval(str(expr), namespace) == expr`` round-trips exactly.  The namespace
is sympy's public names plus a Symbol for each name in the record's
``symbols`` field.  When plain ``str`` does not round-trip, the output of
``EvalSafeStrPrinter`` is tried (rationals as ``S(3)/2``, undefined
functions as ``Function('F')(x)``, constructor form for ``Mul``s that
pairwise evaluation would re-distribute); expressions containing Dummy
variables are accepted when ``dummy_eq`` holds and are listed in the
record's ``dummy_eq_fields`` key.  Fields for which no readable form
round-trips are stored as ``srepr`` instead and named in the record's
``srepr_fields`` key.  The run is resumable: cases whose
result line is already present are skipped, so the script can simply be
re-launched after an interruption.  Files (and cases within each file) are
processed in random order.

Each case runs in a forked child process with a hard timeout, so hangs,
crashes and RecursionErrors in one integral cannot take down the run.

Usage::

    python run_rubi_manualintegrate.py [--timeout 30] [--jobs 7] [--limit N]
    python run_rubi_manualintegrate.py --summarize
"""
from __future__ import annotations

import argparse
import importlib
import json
import multiprocessing
import os
import random
import sys
import time

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MATCHPY_ROOT = os.path.abspath(os.path.join(REPO_ROOT, os.pardir, "matchpy"))
DEFAULT_RESULTS_DIR = os.path.join(REPO_ROOT, "rubi_manualintegrate_results")

# Make sure the sympy from this repository is the one being tested, even if
# another sympy is installed in the active virtualenv.
sys.path.insert(0, REPO_ROOT)

RULE_REPR_MAX = 50000
NUMERIC_TRIES = 8
NUMERIC_NEEDED = 3
NUMERIC_TOL = 1e-8

#: JSON fields that hold serialized SymPy expressions.
EXPR_FIELDS = ("integrand", "variable", "expected_integral",
               "integrand_original", "antiderivative")


def _sympy_namespace(names):
    """Namespace for eval-ing serialized expressions: everything public
    from sympy, plus a plain Symbol for each given name (shadowing any
    sympy names)."""
    import sympy
    ns = {k: v for k, v in vars(sympy).items() if not k.startswith("_")}
    ns.update({name: sympy.Symbol(name) for name in names})
    return ns


def _make_eval_safe_printer():
    """A StrPrinter whose output evals back to the original expression:

    * rationals print as ``S(3)/2`` (``sympy_integers`` setting), so
      Python does not turn ``3/2`` into a float;
    * ``Mul``s whose leading Number coefficient sits next to an ``Add``
      factor print in constructor form ``Mul(S(2), 1/d, c + d*x)``:
      the pairwise evaluation of ``S(2)*(c + d*x)/d`` would re-distribute
      the coefficient over the Add, changing the structure;
    * undefined functions print as ``Function('F')(x)``, so nothing but
      sympy itself is needed in the namespace;
    * ``Dummy`` variables print by bare name (reconstruction can then only
      be checked with ``dummy_eq``).
    """
    from sympy import Add
    from sympy.core.function import AppliedUndef
    from sympy.printing.str import StrPrinter

    class EvalSafeStrPrinter(StrPrinter):
        def _print_Mul(self, expr):
            # A Number coefficient printed next to an Add — directly, or
            # via the coefficient's denominator meeting a 1/(...) factor
            # in the printed denominator — re-distributes on pairwise
            # eval: S(2)*(x + y) becomes 2*x + 2*y.  Constructor form
            # rebuilds the Mul from all args at once, which is stable.
            args = expr.args
            if args and args[0].is_Number and any(
                    a.is_Add or (a.is_Pow and a.base.is_Add
                                 and a.exp.is_negative is True)
                    for a in args[1:]):
                return "Mul(%s)" % ", ".join(self._print(a) for a in args)
            return super()._print_Mul(expr)

        def _print_Function(self, expr):
            if isinstance(expr, AppliedUndef):
                return "Function(%r)(%s)" % (
                    expr.func.__name__,
                    ", ".join(self._print(a) for a in expr.args))
            return super()._print_Function(expr)

        def _print_Dummy(self, expr):
            return expr.name

    return EvalSafeStrPrinter({"sympy_integers": True})


def serialize_expr(expr):
    """Serialize *expr* as a readable string when it round-trips via
    ``eval`` in a sympy namespace with the expression's symbols declared.
    Plain ``str(expr)`` is tried first, then the eval-safe printer; the
    round trip must reproduce the object exactly (identical ``srepr``),
    except for expressions containing ``Dummy`` variables, which are
    accepted when ``dummy_eq`` holds.  ``srepr`` is the last resort.

    Returns ``(string, kind, symbol_names)`` with *kind* one of
    ``"str"`` (exact round trip), ``"dummy_eq"`` (round trip up to dummy
    renaming) or ``"srepr"``.
    """
    from sympy import srepr, Basic, Symbol, Dummy

    # Bound (e.g. Lambda) variables need declaring too, hence atoms() and
    # not free_symbols.
    atoms = expr.atoms(Symbol)
    has_dummy = any(isinstance(s, Dummy) for s in atoms)
    names = sorted({s.name for s in atoms})
    ns = _sympy_namespace(names)
    reference = srepr(expr)

    candidates = []
    if not has_dummy:
        # str() prints Dummy('w') as '_w'; only the eval-safe printer
        # handles dummies.
        candidates.append(str(expr))
    candidates.append(_make_eval_safe_printer().doprint(expr))
    for text in candidates:
        try:
            parsed = eval(text, {"__builtins__": {}}, dict(ns))
        except Exception:
            continue
        if not isinstance(parsed, Basic):
            continue
        if parsed == expr and srepr(parsed) == reference:
            return text, "str", names
        if has_dummy and parsed.dummy_eq(expr):
            return text, "dummy_eq", names
    return reference, "srepr", names


def build_translation_map():
    """Map Mathematica function names (kept as undefined ``Function('...')``
    by the test-suite generator) to real SymPy functions."""
    from sympy import (LambertW, polylog, gamma, uppergamma, expint, Ei,
                       fresnels, fresnelc, erf, erfc, erfi, Si, Ci, Shi, Chi,
                       zeta, li, loggamma, factorial, digamma, polygamma)

    def _gamma(args):
        if len(args) == 1:
            return gamma(args[0])
        if len(args) == 2:
            return uppergamma(args[0], args[1])
        return uppergamma(args[0], args[1]) - uppergamma(args[0], args[2])

    def _polygamma(args):
        if len(args) == 1:
            return digamma(args[0])
        return polygamma(args[0], args[1])

    return {
        "ProductLog": lambda a: LambertW(*a),
        "PolyLog": lambda a: polylog(*a),
        "Gamma": _gamma,
        "ExpIntegralE": lambda a: expint(*a),
        "ExpIntegralEi": lambda a: Ei(*a),
        "FresnelS": lambda a: fresnels(*a),
        "FresnelC": lambda a: fresnelc(*a),
        "Erf": lambda a: erf(*a),
        "Erfc": lambda a: erfc(*a),
        "Erfi": lambda a: erfi(*a),
        "SinIntegral": lambda a: Si(*a),
        "CosIntegral": lambda a: Ci(*a),
        "SinhIntegral": lambda a: Shi(*a),
        "CoshIntegral": lambda a: Chi(*a),
        "PolyGamma": _polygamma,
        "Zeta": lambda a: zeta(*a),
        "LogIntegral": lambda a: li(*a),
        "LogGamma": lambda a: loggamma(*a),
        "Factorial": lambda a: factorial(*a),
    }


def translate_mathematica_functions(expr):
    from sympy.core.function import AppliedUndef

    mapping = build_translation_map()

    def query(e):
        return isinstance(e, AppliedUndef) and e.func.__name__ in mapping

    def repl(e):
        return mapping[e.func.__name__](e.args)

    return expr.replace(query, repl)


def verify_antiderivative(antiderivative, integrand, variable):
    """Check d/dx(antiderivative) - integrand == 0.

    Returns (verified, method) with verified in (True, False, None); None
    means the check was inconclusive (e.g. numeric evaluation failed).
    """
    from sympy import diff, simplify, Rational, nan, oo, zoo

    dd = diff(antiderivative, variable) - integrand
    if dd == 0:
        return True, "syntactic"
    simplified = simplify(dd)
    if simplified == 0:
        return True, "simplify"

    # Numeric fallback: substitute random values for all free symbols and
    # check that the difference evaluates to (approximately) zero.
    free = sorted(simplified.free_symbols, key=lambda s: s.name)
    rng = random.Random(1234)
    successes = 0
    for _ in range(NUMERIC_TRIES):
        subs = {s: Rational(rng.randint(11, 299), 100) for s in free}
        try:
            val = simplified.subs(subs).n(15, chop=True)
        except Exception:
            continue
        if not val.is_number or val.has(nan, oo, zoo) or val.is_finite is False:
            continue
        if abs(val) > NUMERIC_TOL:
            return False, "numeric"
        successes += 1
        if successes >= NUMERIC_NEEDED:
            return True, "numeric"
    return None, "numeric-inconclusive"


def child_worker(conn, integrand, variable):
    """Run one test case; stream partial results through *conn* so that the
    parent still has the ``integral_steps`` info if verification times out."""
    sys.setrecursionlimit(20000)
    from sympy import Integral
    from sympy.integrals.manualintegrate import integral_steps, manualintegrate

    rec = {}

    def send(stage):
        rec["stage"] = stage
        conn.send(dict(rec))

    try:
        rule = integral_steps(integrand, variable)
        rule_repr = repr(rule)
        if len(rule_repr) > RULE_REPR_MAX:
            rule_repr = rule_repr[:RULE_REPR_MAX]
            rec["rule_truncated"] = True
        rec["rule"] = rule_repr
        rec["rule_class"] = type(rule).__name__
        rec["dont_know"] = rule.contains_dont_know()
    except Exception as e:
        rec["steps_error"] = "%s: %s" % (type(e).__name__, e)
    send("steps")

    antiderivative = None
    if "steps_error" not in rec:
        try:
            antiderivative = manualintegrate(integrand, variable)
            text, kind, names = serialize_expr(antiderivative)
            rec["antiderivative"] = text
            rec["antiderivative_kind"] = kind
            rec["antiderivative_symbols"] = names
            rec["has_unevaluated_integral"] = bool(antiderivative.has(Integral))
        except Exception as e:
            rec["integrate_error"] = "%s: %s" % (type(e).__name__, e)
        send("integrated")

    if antiderivative is not None and not rec.get("has_unevaluated_integral"):
        try:
            verified, method = verify_antiderivative(
                antiderivative, integrand, variable)
            rec["verified"] = verified
            rec["verify_method"] = method
        except Exception as e:
            rec["verified"] = None
            rec["verify_method"] = "error: %s: %s" % (type(e).__name__, e)
    send("final")
    conn.close()


def patch_applied_undef_list_args():
    """The generated test suite applies undefined functions to plain
    Python lists, e.g. ``Function('HypergeometricPFQ')([1, 1, 1], ...)``,
    which crashes AppliedUndef construction (lists are unhashable and not
    Basic).  Convert such arguments to ``Tuple`` so the test modules can
    at least be imported."""
    from sympy import Tuple
    from sympy.core.function import AppliedUndef

    if getattr(AppliedUndef, "_list_args_patched", False):
        return
    orig = AppliedUndef.__new__

    def __new__(cls, *args, **options):
        args = tuple(Tuple(*a) if isinstance(a, (list, tuple)) else a
                     for a in args)
        return orig(cls, *args, **options)

    AppliedUndef.__new__ = __new__
    AppliedUndef._list_args_patched = True


def discover_test_modules(matchpy_root):
    suite_dir = os.path.join(matchpy_root, "rubi_rules", "rubi_test_suite")
    modules = []
    for dirpath, dirnames, filenames in os.walk(suite_dir):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in sorted(filenames):
            if fn.startswith("t_") and fn.endswith(".py"):
                rel = os.path.relpath(os.path.join(dirpath, fn), suite_dir)
                mod = "rubi_rules.rubi_test_suite." + \
                    rel[:-3].replace(os.sep, ".")
                modules.append((rel, mod))
    return modules


def load_done_indices(out_path):
    done = set()
    if not os.path.exists(out_path):
        return done
    with open(out_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                done.add(json.loads(line)["index"])
            except (ValueError, KeyError):
                continue
    return done


def reserialize_results(results_dir):
    """Convert existing result files from all-srepr serialization to the
    str-with-roundtrip-check scheme, adding the ``symbols`` (and, when
    needed, ``srepr_fields``) keys.  Already-converted lines (those that
    have a ``symbols`` key) are left untouched.  Files are rewritten
    atomically."""
    from sympy import sympify

    for dirpath, _dirnames, filenames in os.walk(results_dir):
        for fn in sorted(filenames):
            if not fn.endswith(".jsonl"):
                continue
            path = os.path.join(dirpath, fn)
            out_lines = []
            converted = 0
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    if ("symbols" in rec and "functions" not in rec
                            and not rec.get("srepr_fields")):
                        out_lines.append(json.dumps(rec))
                        continue
                    for key in ("symbols", "functions", "srepr_fields",
                                "dummy_eq_fields"):
                        rec.pop(key, None)
                    symbol_names = set()
                    srepr_fields = []
                    dummy_eq_fields = []
                    new_rec = {}
                    for key, value in rec.items():
                        if key not in EXPR_FIELDS:
                            new_rec[key] = value
                            continue
                        # srepr of hyper() is not sympifiable (TupleArg);
                        # patch it so those old records can be recovered.
                        source = value.replace("TupleArg(", "Tuple(")
                        try:
                            expr = sympify(source)
                        except Exception:
                            new_rec[key] = value
                            srepr_fields.append(key)
                            continue
                        text, kind, names = serialize_expr(expr)
                        symbol_names.update(names)
                        if kind == "srepr":
                            srepr_fields.append(key)
                        elif kind == "dummy_eq":
                            dummy_eq_fields.append(key)
                        new_rec[key] = text
                    new_rec["symbols"] = sorted(symbol_names)
                    if dummy_eq_fields:
                        new_rec["dummy_eq_fields"] = sorted(dummy_eq_fields)
                    if srepr_fields:
                        new_rec["srepr_fields"] = sorted(srepr_fields)
                    out_lines.append(json.dumps(new_rec))
                    converted += 1
            if converted:
                tmp_path = path + ".tmp"
                with open(tmp_path, "w") as f:
                    f.write("\n".join(out_lines) + "\n")
                os.replace(tmp_path, path)
            print("%s: %d lines converted" % (path, converted), flush=True)


def summarize(results_dir):
    from collections import Counter, defaultdict
    per_dir = defaultdict(Counter)
    for dirpath, _dirnames, filenames in os.walk(results_dir):
        for fn in filenames:
            if not fn.endswith(".jsonl"):
                continue
            top = os.path.relpath(dirpath, results_dir).split(os.sep)[0]
            with open(os.path.join(dirpath, fn)) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except ValueError:
                        continue
                    c = per_dir[top]
                    c["total"] += 1
                    status = rec.get("status")
                    if status != "ok":
                        c[status if status in ("timeout",) else "failed"] += 1
                    if rec.get("dont_know"):
                        c["dont_know"] += 1
                    solved = (status == "ok" and not rec.get("dont_know")
                              and "antiderivative" in rec
                              and not rec.get("has_unevaluated_integral"))
                    if solved:
                        c["solved"] += 1
                        if rec.get("verified") is True:
                            c["verified"] += 1
                        elif rec.get("verified") is False:
                            c["wrong"] += 1
    total = Counter()
    header = ("directory", "total", "solved", "verified", "wrong",
              "dont_know", "timeout", "failed")
    print(("%-40s" + "%10s" * 7) % header)
    for top in sorted(per_dir):
        c = per_dir[top]
        total.update(c)
        print(("%-40s" + "%10d" * 7) % (top, c["total"], c["solved"],
              c["verified"], c["wrong"], c["dont_know"], c["timeout"],
              c["failed"]))
    print(("%-40s" + "%10d" * 7) % ("TOTAL", total["total"], total["solved"],
          total["verified"], total["wrong"], total["dont_know"],
          total["timeout"], total["failed"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matchpy-root", default=DEFAULT_MATCHPY_ROOT)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--timeout", type=float, default=30.0,
                        help="per-case timeout in seconds")
    parser.add_argument("--jobs", type=int,
                        default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--limit", type=int, default=None,
                        help="stop after this many new cases")
    parser.add_argument("--seed", type=int, default=None,
                        help="seed for the file/case shuffling")
    parser.add_argument("--summarize", action="store_true",
                        help="print a summary of existing results and exit")
    parser.add_argument("--reserialize", action="store_true",
                        help="convert existing all-srepr result files to "
                             "the str-based serialization and exit")
    args = parser.parse_args()

    if args.summarize:
        summarize(args.results_dir)
        return
    if args.reserialize:
        reserialize_results(args.results_dir)
        return

    sys.path.insert(0, args.matchpy_root)
    import sympy
    patch_applied_undef_list_args()
    print("Testing sympy %s from %s" % (sympy.__version__, sympy.__file__),
          flush=True)

    rng = random.Random(args.seed)
    ctx = multiprocessing.get_context("fork")

    modules = discover_test_modules(args.matchpy_root)
    rng.shuffle(modules)

    processed = 0
    t0 = time.time()
    for rel, mod_name in modules:
        out_path = os.path.join(args.results_dir, rel[:-3] + ".jsonl")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        done = load_done_indices(out_path)

        try:
            module = importlib.import_module(mod_name)
        except Exception as e:
            print("!!! failed to import %s: %s: %s" %
                  (rel, type(e).__name__, e), flush=True)
            continue
        cases = module.TEST_CASES
        indices = [i for i in range(len(cases)) if i not in done]
        if not indices:
            continue
        rng.shuffle(indices)
        print("=== %s: %d cases, %d already done ===" %
              (rel, len(cases), len(done)), flush=True)

        pending = list(indices)
        running = []  # list of dicts with proc state
        out = open(out_path, "a")
        try:
            while pending or running:
                while pending and len(running) < args.jobs and (
                        args.limit is None or
                        processed + len(running) < args.limit):
                    i = pending.pop()
                    tc = cases[i]
                    integrand = translate_mathematica_functions(tc.integrand)
                    case = {"integrand": integrand, "variable": tc.variable}
                    parent_conn, child_conn = ctx.Pipe(duplex=False)
                    proc = ctx.Process(
                        target=child_worker,
                        args=(child_conn, integrand, tc.variable))
                    proc.start()
                    child_conn.close()
                    running.append({
                        "proc": proc, "conn": parent_conn, "partial": {},
                        "start": time.time(), "index": i, "tc": tc,
                        "integrand": integrand,
                    })
                if not running:
                    break

                time.sleep(0.05)
                still = []
                for st in running:
                    proc, conn = st["proc"], st["conn"]
                    while conn.poll():
                        try:
                            st["partial"].update(conn.recv())
                        except EOFError:
                            break
                    finished = not proc.is_alive()
                    timed_out = (not finished and
                                 time.time() - st["start"] > args.timeout)
                    if timed_out:
                        proc.kill()
                        st["partial"]["status"] = "timeout"
                    if not (finished or timed_out):
                        still.append(st)
                        continue
                    proc.join()
                    while conn.poll():
                        try:
                            st["partial"].update(conn.recv())
                        except EOFError:
                            break
                    conn.close()
                    rec = st["partial"]
                    if "status" not in rec:
                        if rec.get("stage") == "final":
                            rec["status"] = "ok"
                        elif proc.exitcode not in (0, None):
                            rec["status"] = ("crashed(exitcode=%s)"
                                             % proc.exitcode)
                        else:
                            rec["status"] = "incomplete"
                    rec.pop("stage", None)
                    tc = st["tc"]
                    symbol_names = set()
                    srepr_fields = []
                    dummy_eq_fields = []

                    def ser(field, expr):
                        text, kind, names = serialize_expr(expr)
                        symbol_names.update(names)
                        if kind == "srepr":
                            srepr_fields.append(field)
                        elif kind == "dummy_eq":
                            dummy_eq_fields.append(field)
                        return text

                    record = {
                        "index": st["index"],
                        "source_file": getattr(module, "SOURCE_FILE", rel),
                        "integrand": ser("integrand", st["integrand"]),
                        "variable": ser("variable", tc.variable),
                        "expected_integral": ser("expected_integral",
                                                 tc.integral),
                        "rubi_num_steps": tc.num_steps,
                        "time": round(time.time() - st["start"], 3),
                    }
                    if st["integrand"] != tc.integrand:
                        record["integrand_original"] = ser(
                            "integrand_original", tc.integrand)
                    record.update(rec)
                    symbol_names.update(
                        record.pop("antiderivative_symbols", []))
                    anti_kind = record.pop("antiderivative_kind", None)
                    if anti_kind == "srepr":
                        srepr_fields.append("antiderivative")
                    elif anti_kind == "dummy_eq":
                        dummy_eq_fields.append("antiderivative")
                    record["symbols"] = sorted(symbol_names)
                    if dummy_eq_fields:
                        record["dummy_eq_fields"] = sorted(dummy_eq_fields)
                    if srepr_fields:
                        record["srepr_fields"] = sorted(srepr_fields)
                    out.write(json.dumps(record) + "\n")
                    out.flush()
                    processed += 1
                    solved = (rec["status"] == "ok"
                              and not rec.get("dont_know")
                              and "antiderivative" in rec
                              and not rec.get("has_unevaluated_integral"))
                    print("[%6d] %s #%d: %s%s verified=%s (%.1fs)" % (
                        processed, rel, st["index"],
                        rec.get("rule_class", "-"),
                        " SOLVED" if solved else
                        (" " + rec["status"] if rec["status"] != "ok"
                         else " dont_know"),
                        rec.get("verified"), record["time"]), flush=True)
                    if processed % 100 == 0:
                        rate = processed / (time.time() - t0)
                        print("--- %d cases processed, %.2f cases/s ---" %
                              (processed, rate), flush=True)
                running = still
                if args.limit is not None and processed >= args.limit:
                    break
        finally:
            out.close()
        if args.limit is not None and processed >= args.limit:
            print("Reached --limit of %d cases, stopping." % args.limit,
                  flush=True)
            break

    print("Done: %d new cases in %.1f s." % (processed, time.time() - t0),
          flush=True)


if __name__ == "__main__":
    main()
