# -*-mode: python; flycheck-mode: nil -*-

$XONSH_SHOW_TRACEBACK = True
$RAISE_SUBPROC_ERROR = True

trace on

import os
import sys
import unicodedata
from collections import defaultdict
from collections.abc import Mapping
from getpass import getpass
from contextlib import contextmanager
import json
import glob
import stat
import configparser
import time

import requests
from requests.auth import HTTPBasicAuth

from rever.activity import activity
from rever.conda import run_in_conda_env

cd ..

$ACTIVITIES = [
    # 'version_bump',
    '_version',
    'mailmap_update',
    'test_sympy',
    'source_tarball',
    'wheel',
    'build_docs',
    'copy_release_files',
    'compare_tar_against_git',
    'test_tarball35',
    'test_tarball36',
    'test_tarball37',
    'test_tarball38',
    'test_wheel35',
    'test_wheel36',
    'test_wheel37',
    'test_wheel38',
    'print_authors',
    'sha256',
    # 'tag',
]

version = $VERSION

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
    with run_in_conda_env(['python=3.7', 'mpmath']):
        ./bin/mailmap_update.py

@activity
def test_sympy():
    with run_in_conda_env(['mpmath', 'matplotlib>=2.2', 'numpy', 'scipy', 'aesara',
        'ipython', 'gmpy2', 'symengine', 'libgfortran', 'cython',
        'tensorflow', 'llvmlite', 'wurlitzer', 'autowrap',
        'python-symengine=0.3.*', 'numexpr', 'antlr-python-runtime>=4.7,<4.8',
        'antlr>=4.7,<4.8'], 'sympy-tests'):

        ./setup.py test

@activity(deps={'_version', 'mailmap_update'})
def source_tarball():
    with run_in_conda_env(['mpmath', 'python=3.7'], 'sympy-release'):
        # Assumes this is run in Docker and git is already clean
        ./setup.py sdist --keep-temp


@activity(deps={'_version', 'mailmap_update'})
def wheel():
    with run_in_conda_env(['mpmath', 'python=3.7', 'setuptools', 'pip', 'wheel'], 'sympy-release'):
        # Assumes this is run in Docker and git is already clean
        ./setup.py bdist_wheel --keep-temp

@activity(deps={'_version'})
def build_docs():
    with run_in_conda_env(['sphinx', 'docutils', 'numpy', 'mpmath', 'matplotlib', 'sphinx-math-dollar'],
        envname='sympy-release-docs'):

        cd doc
        make clean
        make html
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


@activity(deps={'source_tarball', 'wheel', 'build_docs'})
def copy_release_files():
    ls dist
    cp dist/* /root/release/

@activity(deps={'source_tarball'})
def test_tarball37():
    test_tarball('3.7')

@activity(deps={'source_tarball'})
def test_tarball38():
    test_tarball('3.8')

@activity(deps={'source_tarball'})
def test_tarball39():
    test_tarball('3.9')

@activity(deps={'wheel'})
def test_wheel37():
    test_wheel('3.7')

@activity(deps={'wheel'})
def test_wheel38():
    test_wheel('3.8')

@activity(deps={'wheel'})
def test_wheel39():
    test_wheel('3.9')

@activity(deps={'source_tarball'})
def compare_tar_against_git():
    """
    Compare the contents of the tarball against git ls-files

    See the bottom of the file for the whitelists.
    """
    release/compare_tar_against_git.py /root/release/@(tarball_format['source']) .

@activity(deps={'source_tarball', 'wheel'})
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

@activity(deps={'source_tarball', 'wheel', 'build_docs'})
def sha256():
    """
    Print the sha256 sums of the release files
    """
    _sha256(print_=True)

def _sha256(print_=True, local=False):
    if local:
        out = $(shasum -a 256 @(release_files()))
    else:
        out = $(shasum -a 256 /root/release/*)
    # Remove the release/ part for printing. Useful for copy-pasting into the
    # release notes.
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out

@activity(deps={'mailmap_update', 'sha256', 'print_authors',
                'source_tarball', 'wheel', 'build_docs',
                'compare_tar_against_git', 'test_tarball35', 'test_tarball36',
                'test_tarball37', 'test_tarball38', 'test_wheel35',
                'test_wheel36', 'test_wheel37', 'test_wheel38', 'test_sympy'})
def release():
    pass

@activity(deps={'_version'})
def GitHub_release():
    # _GitHub_release(version)
    # Run the github_release.py script
    pass

@GitHub_release.undoer
def GitHub_release():
    # Prevent default undo
    pass

@activity(deps={'_version'})
def update_docs():
    _update_docs()


@activity(deps={'_version'})
def update_sympy_org():
    _update_sympy_org()

@activity()
def update_websites():
    _update_docs()
    _update_sympy_org()

# HELPER FUNCTIONS

def test_tarball(py_version):
    """
    Test that the tarball can be unpacked and installed, and that sympy
    imports in the install.
    """
    if py_version not in {'3.7', '3.8', '3.9'}: # TODO: Add win32
        raise ValueError("release must be one of 3.7, 3.8 or 3.9 not %s" % py_version)


    with run_in_conda_env(['python=%s' % py_version], 'test-install-%s' % py_version):
        cp @('/root/release/{source}'.format(**tarball_format)) @("releasetar.tar.gz".format(**tarball_format))
        tar xvf releasetar.tar.gz

        cd @("{source-orig-notar}".format(**tarball_format))
        python setup.py install
        python -c "import sympy; print(sympy.__version__); print('sympy installed successfully')"
        python -m isympy --help
        isympy --help

def test_wheel(py_version):
    """
    Test that the wheel can be installed, and that sympy imports in the install.
    """
    if py_version not in {'3.7', '3.8', '3.9'}: # TODO: Add win32
        raise ValueError("release must be one of 3.7, 3.8 or 3.9 not %s" % py_version)


    with run_in_conda_env(['python=%s' % py_version], 'test-install-%s' % py_version):
        cp @('/root/release/{wheel}'.format(**tarball_format)) @("{wheel}".format(**tarball_format))
        pip install @("{wheel}".format(**tarball_format))

        python -c "import sympy; print(sympy.__version__); print('sympy installed successfully')"
        python -m isympy --help
        isympy --help

def get_tarball_name(file):
    """
    Get the name of a tarball

    file should be one of

    source-orig:       The original name of the source tarball
    source-orig-notar: The name of the untarred directory
    source:            The source tarball (after renaming)
    wheel:             The wheel
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
    elif file == 'wheel':
        name = 'sympy-{version}-py3-none-any.whl'
    else:
        raise ValueError(file + " is not a recognized argument")

    ret = name.format(version=version, type=file,
        extension=doctypename[file])
    return ret

tarball_name_types = {
    'source-orig',
    'source-orig-notar',
    'source',
    'wheel',
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

def release_files():
    """
    Returns the list of local release files
    """
    return glob.glob('release/release-' + $VERSION + '/*')

def red(text):
    return "\033[31m%s\033[0m" % text

def green(text):
    return "\033[32m%s\033[0m" % text

def yellow(text):
    return "\033[33m%s\033[0m" % text

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

    releaseauthors = set($(git --no-pager log @(old_release_tag + '..') --format=%aN).strip().split('\n'))
    priorauthors = set($(git --no-pager log @(old_release_tag) --format=%aN).strip().split('\n'))
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
            # out, and manually make sure `get_previous_version_tag` is correct
            # assert len(parents) == 2, curtag
            curcommit = curtag + "^" # The parent of the tagged commit
        else:
            print(blue("Using {tag} as the tag for the previous "
                "release.".format(tag=curtag)))
            return curtag
    sys.exit(red("Could not find the tag for the previous release."))

def get_sympy_short_version():
    """
    Get the short version of SymPy being released, not including any rc tags
    (like 0.7.3)
    """
    parts = version.split('.')
    # Remove rc tags
    # Handle both 1.0.rc1 and 1.1rc1
    if not parts[-1].isdigit():
        if parts[-1][0].isdigit():
            parts[-1] = parts[-1][0]
        else:
            parts.pop(-1)
    return '.'.join(parts)

def get_location(location):
    """
    Read/save a location from the configuration file.
    """
    locations_file = os.path.expanduser('~/.sympy/sympy-locations')
    config = configparser.SafeConfigParser()
    config.read(locations_file)
    the_location = config.has_option("Locations", location) and config.get("Locations", location)
    if not the_location:
        the_location = input("Where is the SymPy {location} directory? ".format(location=location))
        if not config.has_section("Locations"):
            config.add_section("Locations")
        config.set("Locations", location, the_location)
        save = raw_input("Save this to file [yes]? ")
        if save.lower().strip() in ['', 'y', 'yes']:
            print("saving to ", locations_file)
            with open(locations_file, 'w') as f:
                config.write(f)
    else:
        print("Reading {location} location from config".format(location=location))

    return os.path.abspath(os.path.expanduser(the_location))

def _update_docs(docs_location=None):
    """
    Update the docs hosted at docs.sympy.org
    """
    docs_location = docs_location or get_location("docs")

    print("Docs location:", docs_location)

    current_version = version
    previous_version = get_previous_version_tag().lstrip('sympy-')

    release_dir = os.path.abspath(os.path.expanduser(os.path.join(os.path.curdir, 'release')))
    docs_zip = os.path.abspath(os.path.join(release_dir, 'release-' + version,
        get_tarball_name('html')))

    cd @(docs_location)

    # Check that the docs directory is clean
    git diff --exit-code > /dev/null
    git diff --cached --exit-code > /dev/null

    git pull

    print("Unzipping docs into repo")
    unzip @(docs_zip) > /dev/null
    mv @(get_tarball_name('html-nozip')) latest

    print("Generating indexes")
    ./generate_indexes.py

    print("Committing")
    git add -A latest
    git commit -a -m @('Updating docs to {version}'.format(version=current_version))

    print("Pushing")
    git push origin

    cd @(release_dir)
    cd ..

def _update_sympy_org(website_location=None):
    """
    Update sympy.org

    This just means adding an entry to the news section.
    """
    website_location = website_location or get_location("sympy.github.com")

    release_dir = os.path.abspath(os.path.expanduser(os.path.join(os.path.curdir, 'release')))

    cd @(website_location)

    # Check that the website directory is clean
    git diff --exit-code > /dev/null
    git diff --cached --exit-code > /dev/null

    git pull

    release_date = time.gmtime(os.path.getctime(os.path.join(release_dir,
        'release-' + version, tarball_format['source'])))
    release_year = str(release_date.tm_year)
    release_month = str(release_date.tm_mon)
    release_day = str(release_date.tm_mday)

    with open(os.path.join(website_location, "templates", "index.html"), 'r') as f:
        lines = f.read().split('\n')
        # We could try to use some html parser, but this way is easier
        try:
            news = lines.index(r"    <h3>{% trans %}News{% endtrans %}</h3>")
        except ValueError:
            error("index.html format not as expected")
        lines.insert(news + 2,  # There is a <p> after the news line. Put it
            # after that.
            r"""        <span class="date">{{ datetime(""" + release_year + """, """ + release_month + """, """ + release_day + """) }}</span> {% trans v='""" + version + """' %}Version {{ v }} released{% endtrans %} (<a href="https://github.com/sympy/sympy/wiki/Release-Notes-for-""" + version + """">{% trans %}changes{% endtrans %}</a>)<br/>
    </p><p>""")

    with open(os.path.join(website_location, "templates", "index.html"), 'w') as f:
        print("Updating index.html template")
        f.write('\n'.join(lines))

    print("Generating website pages")
    ./generate

    print("Committing")
    git commit -a -m @('Add {version} to the news'.format(version=version))

    print("Pushing")
    git push origin

    cd @(release_dir)
    cd ..
