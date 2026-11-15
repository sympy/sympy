# Making a SymPy release

This process uses the release workflow in `.github/workflows/release.yml`.
Pushes to `master` build packages and documentation and test installation.
Pushing a release tag triggers the actual PyPI publication, after manual
approval in GitHub. Releases are direct from the master branch.

## Things to review around a release

- Make an issue for the release where discussions can happen.
- Announce upcoming release on mailing list and ask if there are blockers.
- Outstanding bugs, regressions, and release milestone items.
- Python versions:
  - used in CI jobs.
  - `requires-python` in `pyproject.toml`
  - Python version check in `sympy/__init__.py`.
  - ruff's target Python version.
  - Python classifiers in pyproject.toml.
- mpmath and python-flint version bounds.
- Prepare release notes, including highlights, deprecations, and breaking changes.

## Publishing one release to PyPI

Same procedure for release candidate, final release, or patch release.
Make sure local master is fully up to date for each of these steps.

### 1. Choose the version and commit

No version files or workflows need to be edited because it is all handled by
git tags and setuptools-scm. Need to choose what the version is though e.g.
1.15.0rc1 or 1.15.0 or 1.15.1.

### 2. New release notes page

Before a final (not pre-release) release create a new release notes page on the
wiki for the next version. Update the links between the release notes pages.

### 3. Ensure AUTHORS is current

Make and merge a PR to update the AUTHORS file before making the release:

```bash
python bin/mailmap_check.py --update-authors
```

Review the names and email mappings, merge duplicates in `.mailmap` and rerun.
CC people who merged PRs without checking this at the time and admonish them.

### 4. Create and push the tag

Replace `<commit>` with the selected commit's SHA:

```bash
git tag -a sympy-1.15.0 <commit> -m "SymPy 1.15.0"
git push upstream sympy-1.15.0
```

The tag push starts the release workflow. Do not move or reuse a published
release tag.

If something goes wrong in the release workflow and nothing was published then
you need to delete the tag from the sympy/sympy repo to start again:
```
git push upstream --delete refs/tags/sympy-1.15.0
git tag -d sympy-1.15.0
```

### 5. Approve publication

Open the tagged release workflow run in GitHub Actions. Check the package,
documentation, and installation-test results before approving deployment to the
`pypi` environment. The publish step does not wait for all checks because it is
gated on manual approval.

The small number of people who can approve publication will get an email
notification when it is time to do this.

Double check that the correct commit was tagged before approving.

### 6. Verify the PyPI release

Confirm that the expected version appears on the SymPy PyPI project and that
both the wheel and sdist are there. Do a test `pip install sympy==1.15.0` and
`sympy.test(parallel=True)` in a fresh venv.

### 7. Check other releases

- docs site update (final releases only)
- GitHub release
- Author's list
- Checksums

### 8. Announce release on mailing list

- Release highlights from release notes?
- Author list
- Checksums

