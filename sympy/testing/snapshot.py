"""
Classes and functions for carrying out snapshot testing in sympy.
Ideas and code borrowed from pre-existing doctest framework.
"""

import os
import doctest as pdoctest

from sympy.testing.runtests import (
    SymPyDocTestRunner,
    SymPyOutputChecker,
    SymPyTestResults,
)


class SymPySnapshotRunner(SymPyDocTestRunner):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.snapshot_mismatches = []

    def report_failure(self, out, test, example, got):
        # Normal output mismatch

        if got is None:
            got = ""

        if got and not got.endswith("\n"):
            got += "\n"

        want = example.want or ""
        if want and not want.endswith("\n"):
            want += "\n"

        self.snapshot_mismatches.append({
            "filename": test.filename,
            "lineno": example.lineno,
            "source": example.source,
            "want": want,
            "got": got,
        })

        self._fakeout.truncate(0)
        self._fakeout.seek(0)

    def report_unexpected_exception(self, out, test, example, exc_info):
        # Exception mismatch

        got = format_exception_for_snapshot(exc_info)

        want = example.want or ""
        if want and not want.endswith("\n"):
            want += "\n"

        self.snapshot_mismatches.append({
            "filename": test.filename,
            "lineno": example.lineno,
            "source": example.source,
            "want": want,
            "got": got,
        })

        self._fakeout.truncate(0)
        self._fakeout.seek(0)


# Monkeypatch name-mangled methods for SnapshotRunner
monkeypatched_snapshot_methods = [
    'report_failure',
    'report_unexpected_exception',
]

for method in monkeypatched_snapshot_methods:
    oldname = '_DocTestRunner__' + method
    newname = '_SymPySnapshotRunner__' + method
    setattr(SymPySnapshotRunner, newname, getattr(SymPySnapshotRunner, method))


def apply_snapshot_updates(filename, mismatches):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for m in reversed(mismatches):
        lineno = m["lineno"] - 1
        source_lines = m["source"].count("\n") + 1

        out_start = lineno + source_lines
        out_end = out_start

        # go through old expected output
        while out_end < len(lines):
            line = lines[out_end]
            if line.lstrip().startswith(">>>") or line.lstrip().startswith("```"):
                break
            out_end += 1

        new_output = m["got"].rstrip("\n").splitlines(keepends=True)
        new_output = [l + "\n" for l in new_output]

        # replace it with new output
        lines[out_start:out_end] = new_output

    with open(filename, "w", encoding="utf-8") as f:
        f.writelines(lines)


def _run_section(name, text, filename, base_globs, optionflags, start_line):
    """
    Run snapshot tests for one markdown section.
    """
    globs = base_globs.copy()
    globs["__name__"] = "__main__"

    runner = SymPySnapshotRunner(verbose=False, optionflags=optionflags)
    runner._checker = SymPyOutputChecker()

    parser = pdoctest.DocTestParser()
    test = parser.get_doctest(text, globs, name, filename, 0)

    runner.run(test)
    for m in runner.snapshot_mismatches:
        m["lineno"] = start_line + m["lineno"]

    mismatches = getattr(runner, "snapshot_mismatches", [])

    return runner.tries, mismatches


def snapshot_testfile(filename, module_relative=True, name=None, package=None,
             globs=None, verbose=None, report=True, optionflags=0,
             extraglobs=None, raise_on_error=False,
             parser=pdoctest.DocTestParser(), encoding=None, update=False):

    if package and not module_relative:
        raise ValueError("Package may only be specified for module-"
                         "relative paths.")

    # Relativize the path
    text, filename = pdoctest._load_testfile(
        filename, package, module_relative, encoding)

    sections = split_markdown_sections(text)

    if not sections:
        print(f"{filename}: no test sections found")
        return SymPyTestResults(0, 0)

    preamble_name, preamble_text, preamble_start = sections[0]
    test_sections = sections[1:]

    # Base globals (from preamble)
    base_globs = {}

    if preamble_text.strip():
        preamble_test = parser.get_doctest(
            preamble_text, base_globs, "<preamble>", filename, 0
        )

        preamble_runner = SymPySnapshotRunner(verbose=False, optionflags=optionflags)
        preamble_runner._checker = SymPyOutputChecker()
        # run the preamble code, but do not forget the globals set
        preamble_runner.run(preamble_test, clear_globs=False)

        if preamble_runner.snapshot_mismatches:
            for m in preamble_runner.snapshot_mismatches:
                print(f"Preamble failure at line {m['lineno']}")
                print(m["source"].rstrip())
                print("--- expected")
                print(m["want"].rstrip())
                print("+++ got")
                print(m["got"].rstrip())

            return SymPyTestResults(
                len(preamble_runner.snapshot_mismatches),
                preamble_runner.tries,
            )

        # record the globals set in the preamble
        base_globs = preamble_test.globs.copy()

    if not test_sections:
        print(f"{filename}: no test sections after preamble")
        return SymPyTestResults(0, 0)

    total_tries = 0
    total_failures = 0
    all_mismatches = []

    print(f"{filename}: {len(test_sections)} sections")

    # Run sections
    for name, text, start in test_sections:
        print(f"  [{name}]")

        tries, mismatches = _run_section(
            name,
            text,
            filename,
            base_globs,
            optionflags,
            start,
        )

        total_tries += tries
        if mismatches:
            total_failures += len(mismatches)
            all_mismatches.extend(mismatches)

            for m in mismatches:
                print(f"    Mismatch at line {m['lineno']}")
                print(m["source"].rstrip())
                print("    --- expected")
                print(m["want"].rstrip())
                print("    +++ got")
                print(m["got"].rstrip())

    if update and all_mismatches:
        apply_snapshot_updates(filename, all_mismatches)
        return SymPyTestResults(0, total_tries)

    return SymPyTestResults(total_failures, total_tries)


def split_markdown_sections(text):
    """
    Split file into sections based on '# ' headings.

    Returns:
        sections: list of (name, text, start_line)
    """
    lines = text.splitlines(keepends=True)

    sections = []

    current_name = None
    current_lines = []
    current_start = 0

    for i, line in enumerate(lines):
        if line.startswith("# "):
            # Save previous
            if current_name is not None:
                sections.append(
                    (current_name, "".join(current_lines), current_start)
                )

            current_name = line[2:].strip()
            current_lines = []

            # + 1 because doctest uses 1 based line numbers
            current_start = i + 1
        else:
            current_lines.append(line)

    if current_name is not None:
        sections.append(
            (current_name, "".join(current_lines), current_start)
        )

    return sections


def format_exception_for_snapshot(exc_info):
    """
    Format exception like doctest does:
    Keep only first + last line.
    """
    import traceback

    lines = traceback.format_exception_only(exc_info[0], exc_info[1])

    if not lines:
        return ""

    last = lines[-1].rstrip()

    return (
        "Traceback (most recent call last):\n"
        "...\n"
        + last + "\n"
    )


def iter_markdown_files(paths):
    for path in paths:
        if os.path.isdir(path):
            for root, _, files in os.walk(path):
                for f in files:
                    if f.endswith(".md"):
                        yield os.path.join(root, f)
        else:
            if path.endswith(".md"):
                yield path


def run_snapshot_tests(options, args):
    all_files = list(iter_markdown_files(args))
    ok = True

    if not all_files:
        print("No markdown files found.")
        return True

    for filename in all_files:
        print(f"Running snapshot tests: {filename}")
        try:
            failures, attempted = snapshot_testfile(
                filename=filename,
                module_relative=False,
                encoding="utf-8",
                optionflags=(
                    pdoctest.ELLIPSIS |
                    pdoctest.NORMALIZE_WHITESPACE |
                    pdoctest.IGNORE_EXCEPTION_DETAIL |
                    pdoctest.REPORT_UDIFF |
                    pdoctest.FAIL_FAST
                ),
                update=options.update
            )
        except Exception as e:
            print(f"Error while testing {filename}: {e}")
            ok = False
            continue

        if failures and not options.update:
            ok = False
        else:
            print(f"Successfully tested {filename}.")

    return ok
