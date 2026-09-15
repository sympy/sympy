#!/usr/bin/env python
"""Generate a Markdown report from the results of run_rubi_manualintegrate.

Reads the JSONL result files produced by ``run_rubi_manualintegrate.py``
and writes a report with overall statistics, per-section and per-file
breakdowns, rule usage, verification outcomes, wrong results, error
classes and timing information.

Usage::

    python generate_rubi_report.py [--results-dir DIR] [--output FILE]
"""
from __future__ import annotations

import argparse
import json
import os
from collections import Counter, defaultdict

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESULTS_DIR = os.path.join(REPO_ROOT, "rubi_manualintegrate_results")


def load_records(results_dir):
    records = []
    for dirpath, _dirnames, filenames in os.walk(results_dir):
        for fn in sorted(filenames):
            if not fn.endswith(".jsonl"):
                continue
            rel = os.path.relpath(os.path.join(dirpath, fn), results_dir)
            with open(os.path.join(dirpath, fn)) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    rec["_file"] = rel
                    rec["_section"] = rel.split(os.sep)[0]
                    records.append(rec)
    return records


def is_solved(rec):
    return (rec.get("status") == "ok" and not rec.get("dont_know")
            and "antiderivative" in rec
            and not rec.get("has_unevaluated_integral"))


def classify(rec):
    """One outcome label per record."""
    if is_solved(rec):
        v = rec.get("verified")
        if v is True:
            return "solved_verified"
        if v is False:
            return "solved_wrong"
        return "solved_unverified"
    status = rec.get("status", "?")
    if status == "timeout":
        return "timeout"
    if status != "ok":
        return "failed"
    if rec.get("steps_error") or rec.get("integrate_error"):
        return "error"
    return "dont_know"


OUTCOMES = ("solved_verified", "solved_unverified", "solved_wrong",
            "dont_know", "timeout", "error", "failed")
OUTCOME_TITLES = {
    "solved_verified": "solved, verified correct",
    "solved_unverified": "solved, verification inconclusive",
    "solved_wrong": "solved, WRONG antiderivative",
    "dont_know": "not solved (DontKnowRule / unevaluated Integral)",
    "timeout": "timeout",
    "error": "exception from integral_steps/manualintegrate",
    "failed": "worker crashed / incomplete",
}


def pct(n, d):
    return "%.1f%%" % (100.0 * n / d) if d else "-"


def outcome_table(rows, key_title, totals=True):
    """rows: list of (name, Counter-by-outcome).  Returns markdown lines."""
    lines = ["| %s | cases | solved | verified | wrong | dont_know | "
             "timeout | error | solved %% |" % key_title,
             "|---|---|---|---|---|---|---|---|---|"]
    grand = Counter()
    for name, c in rows:
        total = sum(c.values())
        grand.update(c)
        solved = (c["solved_verified"] + c["solved_unverified"]
                  + c["solved_wrong"])
        lines.append(
            "| %s | %d | %d | %d | %d | %d | %d | %d | %s |" % (
                name, total, solved, c["solved_verified"],
                c["solved_wrong"], c["dont_know"], c["timeout"],
                c["error"] + c["failed"], pct(solved, total)))
    if totals and len(rows) > 1:
        total = sum(grand.values())
        solved = (grand["solved_verified"] + grand["solved_unverified"]
                  + grand["solved_wrong"])
        lines.append(
            "| **TOTAL** | **%d** | **%d** | **%d** | **%d** | **%d** | "
            "**%d** | **%d** | **%s** |" % (
                total, solved, grand["solved_verified"],
                grand["solved_wrong"], grand["dont_know"], grand["timeout"],
                grand["error"] + grand["failed"], pct(solved, total)))
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--output", default=None,
                        help="output file (default: REPORT.md inside the "
                             "results directory)")
    args = parser.parse_args()
    output = args.output or os.path.join(args.results_dir, "REPORT.md")

    records = load_records(args.results_dir)
    for rec in records:
        rec["_outcome"] = classify(rec)

    total = len(records)
    by_outcome = Counter(r["_outcome"] for r in records)
    solved_recs = [r for r in records if r["_outcome"].startswith("solved")]

    md = []
    md.append("# manualintegrate on the Rubi integration test suite\n")
    md.append("Results of running `sympy.integrals.manualintegrate` over "
              "the %d integrals of the Rubi (MathematicaSyntaxTestSuite) "
              "test suite, with a 30 s per-case timeout and correctness "
              "verified by differentiating the returned antiderivative "
              "(`simplify`, then numeric spot checks).\n" % total)

    # --- overall ---
    md.append("## Overall outcome\n")
    md.append("| outcome | cases | share |")
    md.append("|---|---|---|")
    for key in OUTCOMES:
        if by_outcome[key]:
            md.append("| %s | %d | %s |" % (
                OUTCOME_TITLES[key], by_outcome[key], pct(by_outcome[key],
                                                          total)))
    solved = len(solved_recs)
    md.append("| **any antiderivative found** | **%d** | **%s** |"
              % (solved, pct(solved, total)))
    md.append("")
    if solved:
        ok = by_outcome["solved_verified"]
        md.append("Of the %d solved integrals, %d (%s) were verified "
                  "correct, %d could not be verified either way, and %d "
                  "were wrong.\n"
                  % (solved, ok, pct(ok, solved),
                     by_outcome["solved_unverified"],
                     by_outcome["solved_wrong"]))

    # --- per section ---
    md.append("## By test-suite section\n")
    per_section = defaultdict(Counter)
    for r in records:
        per_section[r["_section"]][r["_outcome"]] += 1
    md.extend(outcome_table(sorted(per_section.items()), "section"))
    md.append("")

    # --- by Rubi difficulty ---
    md.append("## Solved rate by Rubi step count\n")
    md.append("Number of integration steps Rubi itself needed, as a "
              "difficulty proxy.\n")
    buckets = ((-10**9, 0), (1, 1), (2, 2), (3, 4), (5, 7), (8, 12),
               (13, 10**9))
    per_bucket = defaultdict(Counter)
    for r in records:
        n = r.get("rubi_num_steps") or 0
        for lo, hi in buckets:
            if lo <= n <= hi:
                label = ("<=0" if lo < -10**8 else
                         "%d" % lo if lo == hi else
                         "%d+" % lo if hi > 10**8 else "%d-%d" % (lo, hi))
                per_bucket[(lo, label)][r["_outcome"]] += 1
                break
    md.extend(outcome_table(
        [(label, c) for (_lo, label), c in sorted(per_bucket.items())],
        "Rubi steps"))
    md.append("")

    # --- rules used in solved cases ---
    md.append("## Top-level rule of solved integrals\n")
    md.append("| rule | solved | verified | wrong |")
    md.append("|---|---|---|---|")
    rule_stats = defaultdict(Counter)
    for r in solved_recs:
        rule_stats[r.get("rule_class", "?")][r["_outcome"]] += 1
    for rule, c in sorted(rule_stats.items(),
                          key=lambda kv: -sum(kv[1].values())):
        md.append("| %s | %d | %d | %d |" % (
            rule, sum(c.values()), c["solved_verified"], c["solved_wrong"]))
    md.append("")

    # --- verification methods ---
    md.append("## Verification method of correct results\n")
    md.append("| method | cases |")
    md.append("|---|---|")
    vm = Counter(r.get("verify_method", "?") for r in records
                 if r.get("verified") is True)
    for method, n in vm.most_common():
        md.append("| %s | %d |" % (method, n))
    md.append("")

    # --- wrong results ---
    wrong = [r for r in records if r["_outcome"] == "solved_wrong"]
    md.append("## Wrong antiderivatives (%d)\n" % len(wrong))
    for r in wrong:
        md.append("* `%s`  \n  source: %s #%d, rule `%s`" % (
            r["integrand"], r["source_file"], r["index"],
            r.get("rule_class")))
    md.append("")

    # --- errors ---
    errors = Counter()
    for r in records:
        for key in ("steps_error", "integrate_error"):
            if r.get(key):
                errors[r[key].split(":")[0]] += 1
    if errors:
        md.append("## Exceptions raised\n")
        md.append("| exception | cases |")
        md.append("|---|---|")
        for exc, n in errors.most_common():
            md.append("| %s | %d |" % (exc, n))
        md.append("")

    # --- timing ---
    times = sorted(r.get("time", 0.0) for r in records)
    if times:
        md.append("## Timing\n")
        n = len(times)
        md.append("Total CPU-case time %.1f h; per case: mean %.2f s, "
                  "median %.2f s, p95 %.2f s, max %.1f s.\n" % (
                      sum(times) / 3600.0, sum(times) / n,
                      times[n // 2], times[int(n * 0.95)], times[-1]))
        slowest = sorted((r for r in solved_recs),
                         key=lambda r: -r.get("time", 0))[:5]
        md.append("Slowest solved cases:\n")
        for r in slowest:
            md.append("* %.1f s — `%s` (%s #%d)" % (
                r.get("time", 0), r["integrand"][:90], r["source_file"],
                r["index"]))
        md.append("")

    # --- per-file appendix ---
    md.append("## Appendix: per-file breakdown\n")
    per_file = defaultdict(Counter)
    for r in records:
        per_file[r["_file"]][r["_outcome"]] += 1
    md.extend(outcome_table(sorted(per_file.items()), "file"))
    md.append("")

    with open(output, "w") as f:
        f.write("\n".join(md))
    print("Report on %d records written to %s" % (total, output))


if __name__ == "__main__":
    main()
