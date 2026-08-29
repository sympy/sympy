Title: docs(printing/LLVM): add docstrings for llvmjitcode module

Short summary
- Add docstrings for `LLVMJitCallbackPrinter`, `LLVMJitCode`, `LLVMJitCodeCallback`, and `CodeSignature` in `sympy/printing/llvmjitcode.py` so Sphinx exposes the API via the existing `automodule` entry in `doc/src/modules/printing.rst`.

Rationale
- Issue #28470 reported that the LLVM JIT printing functionality (e.g., `llvm_callable`) was not present in the online docs. The module already had a docstring for `llvm_callable` and the printing docs already include `.. automodule:: sympy.printing.llvmjitcode :members:`. This PR ensures the main public classes also have docstrings so Sphinx will render them.

Testing and CI
- Documentation-only change; no behavior changes. CI should pass; maintainers will build docs on merge.

Release notes (put this exact block into the PR description to satisfy the SymPy release-notes bot)
<!-- BEGIN RELEASE NOTES -->
printing
- Documentation: Add docstrings for LLVM JIT printer classes in sympy.printing.llvmjitcode so Sphinx exposes `llvm_callable` and the LLVM JIT printer classes in the Printing module docs. (doc-only)
<!-- END RELEASE NOTES -->

Closes: #28470

---

Notes for reviewers
- This is a doc-only change. The source-level functionality is unchanged.
- I also added `.mailmap` mapping to include the contributor entry for Amit Joiya so this PR passes the AUTHORS/mailmap check in CI.

How to open a PR
- Copy the above title and description into the PR form when creating the PR.
- Or open this compare link (replace user/repo if necessary):
  https://github.com/sympy/sympy/compare/master...Amitjoiya:docs/llvmjit-printer-docs-pr?expand=1

Thanks!