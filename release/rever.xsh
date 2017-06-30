# -*-mode: python; flycheck-mode: nil -*-

$XONSH_SHOW_TRACEBACK = True
$RAISE_SUBPROC_ERROR = True

import os
import sys
import unicodedata
from collections import defaultdict
from collections.abc import Mapping

from rever.activity import activity
from rever.conda import run_in_conda_env

cd ..

$ACTIVITIES = [
    # 'version_bump',
    '_version',
    'mailmap_update',
    'source_tarball',
    'build_docs',
    'copy_release_files',
    'compare_tar_against_git',
    'test_tarball27',
    'test_tarball33',
    'test_tarball34',
    'test_tarball35',
    'test_tarball36',
    'print_authors',
    # 'tag',
]

# Work around https://github.com/ergs/rever/issues/15
@activity
def _version():
    global version
    version = $VERSION

$TAG_PUSH = False

$VERSION_BUMP_PATTERNS = [
    ('sympy/release.py', r'__version__ = ".*"', r'__version__ = "$VERSION"'),
    ]

# ACTIVITIES

@activity
def mailmap_update():
    with run_in_conda_env(['python=3.6', 'mpmath']):
        ./bin/mailmap_update.py

@activity
def source_tarball():
    with run_in_conda_env(['mpmath', 'python=3.6'], 'sympy-release'):
        # Assumes this is run in Docker and git is already clean
        ./setup.py sdist --keep-temp

@activity
def build_docs():
    with run_in_conda_env(['sphinx=1.3.1', 'docutils=0.12', 'numpy', 'mpmath'],
        envname='sympy-release-docs'):

        cd doc
        make clean
        make html-errors
        make man

        cd _build
        mv html @(tarball_format['html-nozip'])
        zip -9lr @(tarball_format['html']) @(tarball_format['html-nozip'])
        cp @(tarball_format['html']) ../../dist/
        cd ..

        make clean
        make latex

        cd _build/latex
        make
        cp @(tarball_format['pdf-orig']) @("../../../dist/{pdf}".format(**tarball_format))
        cd ../../../


@activity
def copy_release_files():
    ls dist
    cp dist/* /home/release/

@activity
def test_tarball27():
    test_tarball('2.7')

@activity
def test_tarball33():
    test_tarball('3.3')

@activity
def test_tarball34():
    test_tarball('3.4')

@activity
def test_tarball35():
    test_tarball('3.5')

@activity
def test_tarball36():
    test_tarball('3.6')

@activity
def compare_tar_against_git():
    """
    Compare the contents of the tarball against git ls-files

    See the bottom of the file for the whitelists.
    """
    git_lsfiles = set([i.strip() for i in $(git ls-files).strip().split("\n")])
    tar_output_orig = set(show_files('source', print_=False).split("\n"))
    tar_output = set()
    for file in tar_output_orig:
        # The tar files are like sympy-0.7.3/sympy/__init__.py, and the git
        # files are like sympy/__init__.py.
        split_path = full_path_split(file)
        if split_path[-1]:
            # Exclude directories, as git ls-files does not include them
            tar_output.add(os.path.join(*split_path[1:]))
    # print tar_output
    # print git_lsfiles
    fail = False
    print()
    print(blue("Files in the tarball from git that should not be there:"))
    print()
    for line in sorted(tar_output.intersection(git_whitelist)):
        fail = True
        print(line)
    print()
    print(blue("Files in git but not in the tarball:"))
    print()
    for line in sorted(git_lsfiles - tar_output - git_whitelist):
        fail = True
        print(line)
    print()
    print(blue("Files in the tarball but not in git:"))
    print()
    for line in sorted(tar_output - git_lsfiles - tarball_whitelist):
        fail = True
        print(line)
    print()

    if fail:
        sys.exit(red("Non-whitelisted files found or not found in the tarball"))

@activity
def print_authors():
    """
    Print authors text to put at the bottom of the release notes
    """
    authors, authorcount, newauthorcount = get_authors()

    print(blue("Here are the authors to put at the bottom of the release notes."))
    print()
    print("""## Authors

The following people contributed at least one patch to this release (names are
given in alphabetical order by last name). A total of {authorcount} people
contributed to this release. People with a * by their names contributed a
patch for the first time for this release; {newauthorcount} people contributed
for the first time for this release.

Thanks to everyone who contributed to this release!
""".format(authorcount=authorcount, newauthorcount=newauthorcount))

    for name in authors:
        print("- " + name)
    print()

# HELPER FUNCTIONS

def test_tarball(py_version):
    """
    Test that the tarball can be unpacked and installed, and that sympy
    imports in the install.
    """
    if py_version not in {'2.7', '3.3', '3.4', '3.5', '3.6'}: # TODO: Add win32
        raise ValueError("release must be one of 2.7, 3.3, 3.4, 3.5, or 3.6 not %s" % py_version)


    with run_in_conda_env(['python=%s' % py_version], 'test-install-%s' % py_version):
        cp @('/home/release/{source}'.format(**tarball_format)) @("releasetar.tar".format(**tarball_format))
        tar xvf releasetar.tar

        cd @("/home/{source-orig-notar}".format(**tarball_format))
        python setup.py install
        python -c "import sympy; print(sympy.__version__); print('sympy installed successfully')"


def get_tarball_name(file):
    """
    Get the name of a tarball

    file should be one of

    source-orig:       The original name of the source tarball
    source-orig-notar: The name of the untarred directory
    source:            The source tarball (after renaming)
    html:              The name of the html zip
    html-nozip:        The name of the html, without ".zip"
    pdf-orig:          The original name of the pdf file
    pdf:               The name of the pdf file (after renaming)
    """
    doctypename = defaultdict(str, {'html': 'zip', 'pdf': 'pdf'})

    if file in {'source-orig', 'source'}:
        name = 'sympy-{version}.tar.gz'
    elif file == 'source-orig-notar':
        name = "sympy-{version}"
    elif file in {'html', 'pdf', 'html-nozip'}:
        name = "sympy-docs-{type}-{version}"
        if file == 'html-nozip':
            # zip files keep the name of the original zipped directory. See
            # https://github.com/sympy/sympy/issues/7087.
            file = 'html'
        else:
            name += ".{extension}"
    elif file == 'pdf-orig':
        name = "sympy-{version}.pdf"
    else:
        raise ValueError(file + " is not a recognized argument")

    ret = name.format(version=version, type=file,
        extension=doctypename[file])
    return ret

tarball_name_types = {
    'source-orig',
    'source-orig-notar',
    'source',
    'html',
    'html-nozip',
    'pdf-orig',
    'pdf',
    }

# Have to make this lazy so that version can be defined.
class _tarball_format(Mapping):
    def __getitem__(self, name):
        return get_tarball_name(name)

    def __iter__(self):
        return iter(tarball_name_types)

    def __len__(self):
        return len(tarball_name_types)

tarball_format = _tarball_format()

def show_files(file, print_=True):
    """
    Show the contents of a tarball.

    The current options for file are

    source: The source tarball
    html: The html docs zip

    """
    # TODO: Test the unarchived name. See
    # https://github.com/sympy/sympy/issues/7087.
    if file == 'source':
        ret = $(tar tf @("/home/release/{source}".format(**tarball_format)))
    elif file == 'html':
        ret = $(unzip -l @("/home/release/{html}".format(**tarball_format)))
    else:
        raise ValueError(file + " is not valid")
    if print_:
        print(ret)
    return ret

def red(text):
    return "\033[31m%s\033[0m" % text

def green(text):
    return "\033[32m%s\033[0m" % text

def yellow(text):
    return "\033[33m%s\033[0m" % text

def blue(text):
    return "\033[34m%s\033[0m" % text

def full_path_split(path):
    """
    Function to do a full split on a path.
    """
    # Based on http://stackoverflow.com/a/13505966/161801
    rest, tail = os.path.split(path)
    if not rest or rest == os.path.sep:
        return (tail,)
    return full_path_split(rest) + (tail,)

def get_authors():
    """
    Get the list of authors since the previous release

    Returns the list in alphabetical order by last name.  Authors who
    contributed for the first time for this release will have a star appended
    to the end of their names.

    Note: it's a good idea to use ./bin/mailmap_update.py (from the base sympy
    directory) to make AUTHORS and .mailmap up-to-date first before using
    this. fab vagrant release does this automatically.
    """
    def lastnamekey(name):
        """
        Sort key to sort by last name

        Note, we decided to sort based on the last name, because that way is
        fair. We used to sort by commit count or line number count, but that
        bumps up people who made lots of maintenance changes like updating
        mpmath or moving some files around.
        """
        # Note, this will do the wrong thing for people who have multi-word
        # last names, but there are also people with middle initials. I don't
        # know of a perfect way to handle everyone. Feel free to fix up the
        # list by hand.

        text = name.strip().split()[-1].lower()
        # Convert things like Čertík to Certik
        return unicodedata.normalize('NFKD', text).encode('ascii', 'ignore')

    old_release_tag = get_previous_version_tag()

    releaseauthors = set($(git --no-pager log @(old_release_tag + '..') --format="%aN").strip().split('\n'))
    priorauthors = set($(git --no-pager log @(old_release_tag) --format="%aN").strip().split('\n'))
    releaseauthors = {name.strip() for name in releaseauthors if name.strip()}
    priorauthors = {name.strip() for name in priorauthors if name.strip()}
    newauthors = releaseauthors - priorauthors
    starred_newauthors = {name + "*" for name in newauthors}
    authors = releaseauthors - newauthors | starred_newauthors
    return (sorted(authors, key=lastnamekey), len(releaseauthors), len(newauthors))

def get_previous_version_tag():
    """
    Get the version of the previous release
    """
    # We try, probably too hard, to portably get the number of the previous
    # release of SymPy. Our strategy is to look at the git tags.  The
    # following assumptions are made about the git tags:

    # - The only tags are for releases
    # - The tags are given the consistent naming:
    #    sympy-major.minor.micro[.rcnumber]
    #    (e.g., sympy-0.7.2 or sympy-0.7.2.rc1)
    # In particular, it goes back in the tag history and finds the most recent
    # tag that doesn't contain the current short version number as a substring.
    shortversion = get_sympy_short_version()
    curcommit = "HEAD"
    while True:
        curtag = $(git describe --abbrev=0 --tags @(curcommit)).strip()
        if shortversion in curtag:
            # If the tagged commit is a merge commit, we cannot be sure
            # that it will go back in the right direction. This almost
            # never happens, so just error
            parents = $(git rev-list --parents -n 1 @(curtag)).strip().split()
            # rev-list prints the current commit and then all its parents
            # If the tagged commit *is* a merge commit, just comment this
            # out, and make sure `fab vagrant get_previous_version_tag` is correct
            assert len(parents) == 2, curtag
            curcommit = curtag + "^" # The parent of the tagged commit
        else:
            print(blue("Using {tag} as the tag for the previous "
                "release.".format(tag=curtag)))
            return curtag
    sys.exit(red("Could not find the tag for the previous release."))

def get_sympy_version(version_cache=[]):
    """
    Get the full version of SymPy being released (like 0.7.3.rc1)
    """
    if version_cache:
        return version_cache[0]
    version = $(python -c "import sympy;print(sympy.__version__, end='')")
    assert '\n' not in version
    assert ' ' not in version
    assert '\t' not in version
    version_cache.append(version)
    return version

def get_sympy_short_version():
    """
    Get the short version of SymPy being released, not including any rc tags
    (like 0.7.3)
    """
    version = get_sympy_version()
    parts = version.split('.')
    # Remove rc tags
    # Handle both 1.0.rc1 and 1.1rc1
    if not parts[-1].isdigit():
        if parts[-1][0].isdigit():
            parts[-1] = parts[-1][0]
        else:
            parts.pop(-1)
    return '.'.join(parts)

## TARBALL WHITELISTS

# If a file does not end up in the tarball that should, add it to setup.py if
# it is Python, or MANIFEST.in if it is not.  (There is a command at the top
# of setup.py to gather all the things that should be there).

# TODO: Also check that this whitelist isn't growing out of date from files
# removed from git.

# Files that are in git that should not be in the tarball
git_whitelist = {
    # Git specific dotfiles
    '.gitattributes',
    '.gitignore',
    '.mailmap',
    # Travis
    '.travis.yml',
    '.ci/durations.json',
    '.ci/generate_durations_log.sh',
    '.ci/parse_durations_log.py',
    # Code of conduct
    'CODE_OF_CONDUCT.md',
    # Pull request template
    'PULL_REQUEST_TEMPLATE.md',
    # Nothing from bin/ should be shipped unless we intend to install it. Most
    # of this stuff is for development anyway. To run the tests from the
    # tarball, use setup.py test, or import sympy and run sympy.test() or
    # sympy.doctest().
    'bin/adapt_paths.py',
    'bin/ask_update.py',
    'bin/authors_update.py',
    'bin/build_doc.sh',
    'bin/coverage_doctest.py',
    'bin/coverage_report.py',
    'bin/deploy_doc.sh',
    'bin/diagnose_imports',
    'bin/doctest',
    'bin/generate_module_list.py',
    'bin/generate_test_list.py',
    'bin/get_sympy.py',
    'bin/mailmap_update.py',
    'bin/py.bench',
    'bin/strip_whitespace',
    'bin/sympy_time.py',
    'bin/sympy_time_cache.py',
    'bin/test',
    'bin/test_import',
    'bin/test_import.py',
    'bin/test_isolated',
    'bin/test_setup.py',
    'bin/test_travis.sh',
    # The notebooks are not ready for shipping yet. They need to be cleaned
    # up, and preferrably doctested.  See also
    # https://github.com/sympy/sympy/issues/6039.
    'examples/advanced/identitysearch_example.ipynb',
    'examples/beginner/plot_advanced.ipynb',
    'examples/beginner/plot_colors.ipynb',
    'examples/beginner/plot_discont.ipynb',
    'examples/beginner/plot_gallery.ipynb',
    'examples/beginner/plot_intro.ipynb',
    'examples/intermediate/limit_examples_advanced.ipynb',
    'examples/intermediate/schwarzschild.ipynb',
    'examples/notebooks/density.ipynb',
    'examples/notebooks/fidelity.ipynb',
    'examples/notebooks/fresnel_integrals.ipynb',
    'examples/notebooks/qubits.ipynb',
    'examples/notebooks/sho1d_example.ipynb',
    'examples/notebooks/spin.ipynb',
    'examples/notebooks/trace.ipynb',
    'examples/notebooks/README.txt',
    # This stuff :)
    'release/.gitignore',
    'release/README.md',
    'release/Vagrantfile',
    'release/fabfile.py',
    'release/Dockerfile',
    'release/Dockerfile-base',
    'release/release.sh',
    'release/rever.xsh',
    # This is just a distribute version of setup.py. Used mainly for setup.py
    # develop, which we don't care about in the release tarball
    'setupegg.py',
    # Example on how to use tox to test Sympy. For development.
    'tox.ini.sample',
    # pytest stuff
    'conftest.py',
    # Encrypted deploy key for deploying dev docs to GitHub
    'github_deploy_key.enc',
    }

# Files that should be in the tarball should not be in git

tarball_whitelist = {
    # Generated by setup.py. Contains metadata for PyPI.
    "PKG-INFO",
    # Generated by setuptools. More metadata.
    'setup.cfg',
    'sympy.egg-info/PKG-INFO',
    'sympy.egg-info/SOURCES.txt',
    'sympy.egg-info/dependency_links.txt',
    'sympy.egg-info/requires.txt',
    'sympy.egg-info/top_level.txt',
    'sympy.egg-info/not-zip-safe',
    }
