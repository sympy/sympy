import doctest as pdoctest
import os
import sys

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
        # Record mismatches
        self.snapshot_mismatches.append({
            "filename": test.filename,
            "lineno": example.lineno,
            "source": example.source,
            "want": example.want,
            "got": got,
        })

        self._fakeout.truncate(0)
        self._fakeout.seek(0)


# Monkeypatch name-mangled methods for SnapshotRunner
monkeypatched_snapshot_methods = [
    'report_failure',
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


# This function is same as sympytestfile, with minor changes adapted to snapshot testing.
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

    # If no name was given, then use the file's name.
    if name is None:
        name = os.path.basename(filename)

    # Assemble the globals.
    if globs is None:
        globs = {}
    else:
        globs = globs.copy()
    if extraglobs is not None:
        globs.update(extraglobs)
    if '__name__' not in globs:
        globs['__name__'] = '__main__'

    if raise_on_error:
        runner = pdoctest.DebugRunner(verbose=verbose, optionflags=optionflags)
    else:
        runner = SymPySnapshotRunner(verbose=verbose, optionflags=optionflags)
        runner._checker = SymPyOutputChecker()

    # Read the file, convert it to a test, and run it.
    test = parser.get_doctest(text, globs, name, filename, 0)
    print(f"{filename}: {len(test.examples)} examples")
    runner.run(test)

    print(f"Tried: {runner.tries}")
    print(f"Snapshot mismatches: {len(runner.snapshot_mismatches)}")

    if report:
        runner.summarize()

    if pdoctest.master is None:
        pdoctest.master = runner
    else:
        pdoctest.master.merge(runner)

    if hasattr(runner, "snapshot_mismatches"):
        for m in runner.snapshot_mismatches:
            print(f"Snapshot mismatch in {m['filename']}:{m['lineno']}")
            print(m["source"].rstrip())
            print("--- expected")
            print(m["want"].rstrip())
            print("+++ got")
            print(m["got"].rstrip())

    mismatches = getattr(runner, "snapshot_mismatches", [])

    if update and mismatches:
        apply_snapshot_updates(filename, mismatches)
        return SymPyTestResults(0, runner.tries)

    return SymPyTestResults(len(mismatches), runner.tries)


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
        sys.exit(0)

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
                    pdoctest.IGNORE_EXCEPTION_DETAIL
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
            print(f"Succesfully tested {filename}.")

    return ok
