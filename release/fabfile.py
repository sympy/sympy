# -*- coding: utf-8 -*-
"""
Fab file for releasing

Please read the README in this directory.

Guide for this file
===================

Vagrant is a tool that gives us a reproducible VM, and fabric is a tool that
we use to run commands on that VM.

Each function in this file should be run as

fab vagrant func

Even those functions that do not use vagrant must be run this way, because of
the vagrant configuration at the bottom of this file.

Any function that should be made avaiable from the command line needs to have
the @task decorator.

Save any files that should be reset between runs somewhere in the repos
directory, so that the remove_userspace() function will clear it.  It's best
to do a complete vagrant destroy before a full release, but that takes a
while, so the remove_userspace() ensures that things are mostly reset for
testing.

Do not enforce any naming conventions on the release branch. By tradition, the
name of the release branch is the same as the version being released (like
0.7.3), but this is not required. Use get_sympy_version() and
get_sympy_short_version() to get the SymPy version (the SymPy __version__
*must* be changed in __init__.py for this to work).
"""
from __future__ import print_function

from collections import defaultdict, OrderedDict

from contextlib import contextmanager

from fabric.api import env, local, run, sudo, cd, hide, task
from fabric.contrib.files import exists
from fabric.colors import blue
from fabric.utils import error

import unicodedata

import os.path

try:
    # https://pypi.python.org/pypi/fabric-virtualenv/
    from fabvenv import virtualenv, make_virtualenv
    # Note, according to fabvenv docs, always use an absolute path with
    # virtualenv().
except ImportError:
    error("fabvenv is required. See https://pypi.python.org/pypi/fabric-virtualenv/")

# Note, it's actually good practice to use absolute paths
# everywhere. Otherwise, you will get surprising results if you call one
# function from another, because your current working directory will be
# whatever it was in the calling function, not ~.  Also, due to what should
# probably be considered a bug, ~ is not treated as an absolute path. You have
# to explicitly write out /home/vagrant/

env.use_ssh_config = True

try:
    # Only works in newer versions of fabric
    env.colorize_errors = True
except AttributeError:
    pass

def full_path_split(path):
    """
    Function to do a full split on a path.
    """
    # Based on http://stackoverflow.com/a/13505966/161801
    rest, tail = os.path.split(path)
    if not rest or rest == os.path.sep:
        return (tail,)
    return full_path_split(rest) + (tail,)

@contextmanager
def use_venv(pyversion):
    """
    Change make_virtualenv to use a given cmd

    pyversion should be '2' or '3'
    """
    pyversion = str(pyversion)
    if pyversion == '2':
        yield
    elif pyversion == '3':
        oldvenv = env.virtualenv
        env.virtualenv = 'virtualenv -p /usr/bin/python3'
        yield
        env.virtualenv = oldvenv
    else:
        raise ValueError("pyversion must be one of '2' or '3', not %s" % pyversion)

@task
def prepare():
    """
    Setup the VM

    This only needs to be run once.  It downloads all the necessary software,
    and a git cache. To reset this, use vagrant destroy and vagrant up.  Note,
    this may take a while to finish, depending on your internet connection
    speed.
    """
    prepare_apt()
    checkout_cache()

@task
def prepare_apt():
    """
    Download software from apt

    Note, on a slower internet connection, this will take a while to finish,
    because it has to download many packages, include latex and all its
    dependencies.
    """
    sudo("apt-get -qq update")
    sudo("apt-get -y install git python3 make python-virtualenv zip python-dev")
    # Needed to build the docs
    sudo("apt-get -y install graphviz inkscape texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra")
    # Our Ubuntu is too old to include Python 3.3
    sudo("apt-get -y install python-software-properties")
    sudo("add-apt-repository -y ppa:fkrull/deadsnakes")
    sudo("apt-get -y update")
    sudo("apt-get -y install python3.3")

@task
def remove_userspace():
    """
    Deletes (!) the SymPy changes. Use with great care.

    This should be run between runs to reset everything.
    """
    run("rm -rf repos")

@task
def checkout_cache():
    """
    Checkout a cache of SymPy

    This should only be run once. The cache is use as a --reference for git
    clone.  This makes deleting and recreating the SymPy a la
    remove_userspace() and gitrepos() and clone very fast.
    """
    run("rm -rf sympy-cache.git")
    run("git clone --bare https://github.com/sympy/sympy.git sympy-cache.git")

@task
def gitrepos(branch=None, fork='sympy'):
    """
    Clone the repo

    fab vagrant prepare (namely, checkout_cache()) must be run first. By
    default, the branch checked out is the same one as the one checked out
    locally. The master branch is not allowed--use a release branch (see the
    README). No naming convention is put on the release branch.

    To test the release, create a branch in your fork, and set the fork
    option.
    """
    with cd("/home/vagrant"):
        if not exists("sympy-cache.git"):
            error("Run fab vagrant prepare first")
    if not branch:
        # Use the current branch (of this git repo, not the one in Vagrant)
        branch = local("git rev-parse --abbrev-ref HEAD", capture=True)
    if branch == "master":
        raise Exception("Cannot release from master")
    run("mkdir -p repos")
    with cd("/home/vagrant/repos"):
        run("git clone --reference ../sympy-cache.git https://github.com/{fork}/sympy.git".format(fork=fork))
        with cd("/home/vagrant/repos/sympy"):
            run("git checkout -t origin/%s" % branch)

@task
def get_sympy_version(version_cache=[]):
    """
    Get the full version of SymPy being released (like 0.7.3.rc1)
    """
    if version_cache:
        return version_cache[0]
    if not exists("/home/vagrant/repos/sympy"):
        gitrepos()
    with cd("/home/vagrant/repos/sympy"):
        version = run('python -c "import sympy;print(sympy.__version__)"')
    assert '\n' not in version
    assert ' ' not in version
    assert '\t' not in version
    version_cache.append(version)
    return version

@task
def get_sympy_short_version():
    """
    Get the short version of SymPy being released, not including any rc tags
    (like 0.7.3)
    """
    version = get_sympy_version()
    return '.'.join(version.split('.')[:3]) # Remove any rc tags

@task
def test_sympy():
    """
    Run the SymPy test suite
    """
    with cd("/home/vagrant/repos/sympy"):
        run("./setup.py test")

@task
def test_tarball(release='2'):
    """
    Test that the tarball can be unpacked and installed, and that sympy
    imports in the install.
    """
    if release not in {'2', '3'}: # TODO: Add win32
        raise ValueError("release must be one of '2', '3', not %s" % release)

    venv = "/home/vagrant/repos/test-{release}-virtualenv".format(release=release)
    tarball_formatter_dict = tarball_formatter()

    with use_venv(release):
        make_virtualenv(venv)
        with virtualenv(venv):
            run("cp /vagrant/release/{source} releasetar.tar".format(**tarball_formatter_dict))
            run("tar xvf releasetar.tar")
            with cd("/home/vagrant/{source-orig-notar}".format(**tarball_formatter_dict)):
                run("python setup.py install")
                run('python -c "import sympy; print(sympy.__version__)"')

@task
def release(branch=None, fork='sympy'):
    """
    Perform all the steps required for the release, except uploading

    In particular, it builds all the release files, and puts them in the
    release/ directory in the same directory as this one.  At the end, it
    prints some things that need to be pasted into various places as part of
    the release.

    To test the release, push a branch to your fork on GitHub and set the fork
    option to your username.
    """
    remove_userspace()
    gitrepos(branch, fork)
    # This has to be run locally because it itself uses fabric. I split it out
    # into a separate script so that it can be used without vagrant.
    local("../bin/mailmap_update.py")
    source_tarball()
    build_docs()
    copy_release_files()
    test_tarball('2')
    test_tarball('3')
    compare_tar_against_git()
    print_authors()
    GitHub_release()

@task
def source_tarball():
    """
    Build the source tarball
    """
    with cd("/home/vagrant/repos/sympy"):
        run("git clean -dfx")
        run("./setup.py clean")
        run("./setup.py sdist")
        run("./setup.py bdist_wininst")
        run("mv dist/{win32-orig} dist/{win32}".format(**tarball_formatter()))

@task
def build_docs():
    """
    Build the html and pdf docs
    """
    with cd("/home/vagrant/repos/sympy"):
        run("mkdir -p dist")
        venv = "/home/vagrant/docs-virtualenv"
        make_virtualenv(venv, dependencies=['sphinx==1.1.3', 'numpy'])
        with virtualenv(venv):
            with cd("/home/vagrant/repos/sympy/doc"):
                run("make clean")
                run("make html-errors")
                with cd("/home/vagrant/repos/sympy/doc/_build"):
                    run("mv html {html-nozip}".format(**tarball_formatter()))
                    run("zip -9lr {html} {html-nozip}".format(**tarball_formatter()))
                    run("cp {html} ../../dist/".format(**tarball_formatter()))
                run("make clean")
                run("make latex")
                with cd("/home/vagrant/repos/sympy/doc/_build/latex"):
                    run("make")
                    run("cp {pdf-orig} ../../../dist/{pdf}".format(**tarball_formatter()))

@task
def copy_release_files():
    """
    Move the release files from the VM to release/ locally
    """
    with cd("/home/vagrant/repos/sympy"):
        run("mkdir -p /vagrant/release")
        run("cp dist/* /vagrant/release/")

@task
def show_files(file, print_=True):
    """
    Show the contents of a tarball.

    The current options for file are

    source: The source tarball
    win: The Python 2 Windows installer (Not yet implemented!)
    html: The html docs zip

    Note, this runs locally, not in vagrant.
    """
    # TODO: Test the unarchived name. See
    # https://code.google.com/p/sympy/issues/detail?id=3988.
    if file == 'source':
        ret = local("tar tf release/{source}".format(**tarball_formatter()), capture=True)
    elif file == 'win':
        # TODO: Windows
        raise NotImplementedError("Windows installers")
    elif file == 'html':
        ret = local("unzip -l release/{html}".format(**tarball_formatter()), capture=True)
    else:
        raise ValueError(file + " is not valid")
    if print_:
        print(ret)
    return ret

# If a file does not end up in the tarball that should, add it to setup.py if
# it is Python, or MANIFEST.in if it is not.  (There is a command at the top
# of setup.py to gather all the things that should be there).

# TODO: Also check that this whitelist isn't growning out of date from files
# removed from git.

# TODO: Address the "why?" comments below.

# Files that are in git that should not be in the tarball
git_whitelist = {
    # Git specific dotfiles
    '.gitattributes',
    '.gitignore',
    '.mailmap',
    # Travis
    '.travis.yml',
    # This is the file you should edit if not enough ends up in the tarball
    'MANIFEST.in',
    # Experimental Cythonization support. Not for production
    'Makefile',
    # Nothing from bin/ should be shipped unless we intend to install it. Most
    # of this stuff is for development anyway. To run the tests from the
    # tarball, use setup.py test, or import sympy and run sympy.test() or
    # sympy.doctest().
    'bin/adapt_paths.py',
    'bin/ask_update.py',
    'bin/coverage_doctest.py',
    'bin/coverage_report.py',
    'bin/doctest',
    'bin/generate_test_list.py',
    'bin/get_sympy.py',
    'bin/py.bench',
    'bin/mailmap_update.py',
    'bin/strip_whitespace',
    'bin/sympy_time.py',
    'bin/sympy_time_cache.py',
    'bin/test',
    'bin/test_import',
    'bin/test_import.py',
    'bin/test_isolated',
    'bin/test_travis.sh',
    # This is also related to Cythonization
    'build.py',
    # The notebooks are not ready for shipping yet. They need to be cleaned
    # up, and preferrably doctested.  See also
    # https://code.google.com/p/sympy/issues/detail?id=2940.
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
    # This stuff :)
    'release/.gitignore',
    'release/README.md',
    'release/Vagrantfile',
    'release/fabfile.py',
    # This is just a distribute version of setup.py. Used mainly for setup.py
    # develop, which we don't care about in the release tarball
    'setupegg.py',
    # We don't ship the benchmarks (why?)
    'sympy/benchmarks/bench_meijerint.py',
    'sympy/benchmarks/bench_symbench.py',
    'sympy/core/benchmarks/bench_arit.py',
    'sympy/core/benchmarks/bench_assumptions.py',
    'sympy/core/benchmarks/bench_basic.py',
    'sympy/core/benchmarks/bench_expand.py',
    'sympy/core/benchmarks/bench_numbers.py',
    'sympy/core/benchmarks/bench_sympify.py',
    'sympy/functions/elementary/benchmarks/bench_exp.py',
    'sympy/functions/special/benchmarks/bench_special.py',
    # More benchmarks
    'sympy/integrals/benchmarks/bench_integrate.py',
    'sympy/integrals/benchmarks/bench_trigintegrate.py',
    'sympy/logic/benchmarks/input/10.cnf',
    'sympy/logic/benchmarks/input/100.cnf',
    'sympy/logic/benchmarks/input/105.cnf',
    'sympy/logic/benchmarks/input/110.cnf',
    'sympy/logic/benchmarks/input/115.cnf',
    'sympy/logic/benchmarks/input/120.cnf',
    'sympy/logic/benchmarks/input/125.cnf',
    'sympy/logic/benchmarks/input/130.cnf',
    'sympy/logic/benchmarks/input/135.cnf',
    'sympy/logic/benchmarks/input/140.cnf',
    'sympy/logic/benchmarks/input/145.cnf',
    'sympy/logic/benchmarks/input/15.cnf',
    'sympy/logic/benchmarks/input/150.cnf',
    'sympy/logic/benchmarks/input/20.cnf',
    'sympy/logic/benchmarks/input/25.cnf',
    'sympy/logic/benchmarks/input/30.cnf',
    'sympy/logic/benchmarks/input/35.cnf',
    'sympy/logic/benchmarks/input/40.cnf',
    'sympy/logic/benchmarks/input/45.cnf',
    'sympy/logic/benchmarks/input/50.cnf',
    'sympy/logic/benchmarks/input/55.cnf',
    'sympy/logic/benchmarks/input/60.cnf',
    'sympy/logic/benchmarks/input/65.cnf',
    'sympy/logic/benchmarks/input/70.cnf',
    'sympy/logic/benchmarks/input/75.cnf',
    'sympy/logic/benchmarks/input/80.cnf',
    'sympy/logic/benchmarks/input/85.cnf',
    'sympy/logic/benchmarks/input/90.cnf',
    'sympy/logic/benchmarks/input/95.cnf',
    'sympy/logic/benchmarks/run-solvers.py',
    'sympy/logic/benchmarks/test-solver.py',
    'sympy/matrices/benchmarks/bench_matrix.py',
    # More benchmarks...
    'sympy/polys/benchmarks/__init__.py',
    'sympy/polys/benchmarks/bench_galoispolys.py',
    'sympy/polys/benchmarks/bench_groebnertools.py',
    'sympy/polys/benchmarks/bench_solvers.py',
    'sympy/series/benchmarks/bench_limit.py',
    'sympy/solvers/benchmarks/bench_solvers.py',
    # Example on how to use tox to test Sympy. For development.
    'tox.ini.sample',
    }

# Files that should be in the tarball should not be in git

tarball_whitelist = {
    "PKG-INFO", # Generated by setup.py. Contains metadata for PyPI.
    }

@task
def compare_tar_against_git():
    """
    Compare the contents of the tarball against git ls-files
    """
    with hide("commands"):
        with cd("/home/vagrant/repos/sympy"):
            git_lsfiles = set([i.strip() for i in run("git ls-files").split("\n")])
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
    print(blue("Files in the tarball from git that should not be there:",
        bold=True))
    print()
    for line in sorted(tar_output.intersection(git_whitelist)):
        fail = True
        print(line)
    print()
    print(blue("Files in git but not in the tarball:", bold=True))
    print()
    for line in sorted(git_lsfiles - tar_output - git_whitelist):
        fail = True
        print(line)
    print()
    print(blue("Files in the tarball but not in git:", bold=True))
    print()
    for line in sorted(tar_output - git_lsfiles - tarball_whitelist):
        fail = True
        print(line)

    if fail:
        error("Non-whitelisted files found or not found in the tarball")

@task
def md5(file='*', print_=True):
    """
    Print the md5 sums of the release files
    """
    out = local("md5sum release/" + file, capture=True)
    # Remove the release/ part for printing. Useful for copy-pasting into the
    # release notes.
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out

descriptions = OrderedDict([
    ('source', "The SymPy source installer.",),
    ('win32', "Python Windows 32-bit installer.",),
    ('html', '''Html documentation for the Python 2 version. This is the same as
the <a href="http://docs.sympy.org/latest/index.html">online documentation</a>.''',),
    ('pdf', '''Pdf version of the <a href="http://docs.sympy.org/latest/index.html"> html documentation</a>.''',),
    ])

@task
def table():
    """
    Make an html table of the downloads.

    This is for pasting into the GitHub releases page. See GitHub_release().
    """
    tarball_formatter_dict = tarball_formatter()
    shortversion = get_sympy_short_version()

    tarball_formatter_dict['version'] = shortversion

    md5s = [i.split('\t') for i in md5(print_=False).split('\n')]
    md5s_dict = {name: md5 for md5, name in md5s}

    table = []

    # http://docs.python.org/2/library/contextlib.html#contextlib.contextmanager. Not
    # recommended as a real way to generate html, but it works better than
    # anything else I've tried.
    @contextmanager
    def tag(name):
        table.append("<%s>" % name)
        yield
        table.append("</%s>" % name)

    with tag('table'):
        with tag('tr'):
            for headname in ["Filename", "Description", "md5"]:
                with tag("th"):
                    table.append(headname)

        for key in descriptions:
            name = get_tarball_name(key)
            with tag('tr'):
                with tag('td'):
                    # code renders better than tt or pre
                    with tag('code'):
                        table.append(name)
                with tag('td'):
                    table.append(descriptions[key].format(**tarball_formatter_dict))
                with tag('td'):
                    table.append(md5s_dict[name])

    out = ' '.join(table)
    return out

@task
def GitHub_release():
    """
    Generate text to put in the GitHub release Markdown box
    """
    shortversion = get_sympy_short_version()
    htmltable = table()
    out = """\
See https://github.com/sympy/sympy/wiki/release-notes-for-{shortversion} for the release notes.

{htmltable}

**Note**: Do not download the `Source code (zip)` or the `Source code (tar.gz)`
files below.
"""
    out = out.format(shortversion=shortversion, htmltable=htmltable)
    print(blue("Here are the release notes to copy into the GitHub release "
        "Markdown form:", bold=True))
    print()
    print(out)
    return out

@task
def get_tarball_name(file):
    """
    Get the name of a tarball

    file should be one of

    source-orig:       The original name of the source tarball
    source-orig-notar: The name of the untarred directory
    source:            The source tarball (after renaming)
    win32-orig:        The original name of the win32 installer
    win32:             The name of the win32 installer (after renaming)
    html:              The name of the html zip
    html-nozip:        The name of the html, without ".zip"
    pdf-orig:          The original name of the pdf file
    pdf:               The name of the pdf file (after renaming)
    """
    version = get_sympy_version()
    doctypename = defaultdict(str, {'html': 'zip', 'pdf': 'pdf'})
    winos = defaultdict(str, {'win32': 'win32', 'win32-orig': 'linux-i686'})

    if file in {'source-orig', 'source'}:
        name = 'sympy-{version}.tar.gz'
    elif file == 'source-orig-notar':
        name = "sympy-{version}"
    elif file in {'win32', 'win32-orig'}:
        name = "sympy-{version}.{wintype}.exe"
    elif file in {'html', 'pdf', 'html-nozip'}:
        name = "sympy-docs-{type}-{version}"
        if not file.endswith('nozip'):
            name += ".{extension}"
    elif file == 'pdf-orig':
        name = "sympy-{version}.pdf"
    else:
        raise ValueError(file + " is not a recognized argument")

    ret = name.format(version=version, type=file,
        extension=doctypename[file], wintype=winos[file])
    return ret

tarball_name_types = {
    'source-orig',
    'source-orig-notar',
    'source',
    'win32-orig',
    'win32',
    'html',
    'html-nozip',
    'pdf-orig',
    'pdf',
    }

# This has to be a function, because you cannot call any function here at
# import time (before the vagrant() function is run).
def tarball_formatter():
    return {name: get_tarball_name(name) for name in tarball_name_types}

@task
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
    with cd("/home/vagrant/repos/sympy"):
        while True:
            curtag = run("git describe --abbrev=0 --tags " +
                curcommit).strip()
            if shortversion in curtag:
                # If the tagged commit is a merge commit, we cannot be sure
                # that it will go back in the right direction. This almost
                # never happens, so just error
                parents = local("git rev-list --parents -n 1 " + curtag,
                    capture=True).strip().split()
                # rev-list prints the current commit and then all its parents
                assert len(parents) == 2, curtag
                curcommit = curtag + "^" # The parent of the tagged commit
            else:
                print(blue("Using {tag} as the tag for the previous "
                    "release.".format(tag=curtag), bold=True))
                return curtag
        error("Could not find the tag for the previous release.")

@task
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

        # Note, you must call unicode() *before* lower, or else it won't
        # lowercase non-ASCII characters like Č -> č
        text = unicode(name.strip().split()[-1], encoding='utf-8').lower()
        # Convert things like Čertík to Certik
        return unicodedata.normalize('NFKD', text).encode('ascii', 'ignore')

    old_release_tag = get_previous_version_tag()
    with cd("/home/vagrant/repos/sympy"), hide('commands'):
        releaseauthors = set(run('git --no-pager log {tag}.. --format="%aN"'.format(tag=old_release_tag)).strip().split('\n'))
        priorauthors = set(run('git --no-pager log {tag} --format="%aN"'.format(tag=old_release_tag)).strip().split('\n'))
        releaseauthors = {name.strip() for name in releaseauthors if name.strip()}
        priorauthors = {name.strip() for name in priorauthors if name.strip()}
        newauthors = releaseauthors - priorauthors
        starred_newauthors = {name + "*" for name in newauthors}
        authors = releaseauthors - newauthors | starred_newauthors
        return (sorted(authors, key=lastnamekey), len(releaseauthors), len(newauthors))

@task
def print_authors():
    """
    Print authors text to put at the bottom of the release notes
    """
    authors, authorcount, newauthorcount = get_authors()

    print(blue("Here are the authors to put at the bottom of the release "
        "notes.", bold=True))
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

# ------------------------------------------------
# PyPI

@task
def upload():
    """
    Upload the files everywhere

    For now, it is just PyPI, because GitHub doesn't seem to have an API.

    """
    distutils_check()
    #pypi_register()
    pypi_upload()

@task
def distutils_check():
    """
    Runs setup.py check
    """
    with cd("/home/vagrant/repos/sympy"):
        run("python setup.py check")
        run("python3 setup.py check")

@task
def pypi_register():
    """
    Register a release with PyPI

    This should only be done for the final release. You need PyPI
    authentication to do this.
    """
    with cd("/home/vagrant/repos/sympy"):
        run("python setup.py register")

@task
def pypi_upload():
    """
    Upload files to PyPI
    """
    with cd("/home/vagrant/repos/sympy"):
        # XXX: Doesn't actually work yet
        run("python setupegg.py upload")

# ------------------------------------------------
# Vagrant related configuration

@task
def vagrant():
    """
    Run commands using vagrant
    """
    vc = get_vagrant_config()
    # change from the default user to 'vagrant'
    env.user = vc['User']
    # connect to the port-forwarded ssh
    env.hosts = ['%s:%s' % (vc['HostName'], vc['Port'])]
    # use vagrant ssh key
    env.key_filename = vc['IdentityFile'].strip('"')
    # Forward the agent if specified:
    env.forward_agent = vc.get('ForwardAgent', 'no') == 'yes'

def get_vagrant_config():
    """
    Parses vagrant configuration and returns it as dict of ssh parameters
    and their values
    """
    result = local('vagrant ssh-config', capture=True)
    conf = {}
    for line in iter(result.splitlines()):
        parts = line.split()
        conf[parts[0]] = ' '.join(parts[1:])
    return conf

@task
def restart_network():
    """
    Do this if the VM won't connect to the internet.
    """
    run("sudo /etc/init.d/networking restart")

# ---------------------------------------
# Just a simple testing command:

@task
def uname():
    """
    Get the uname in Vagrant. Useful for testing that Vagrant works.
    """
    run('uname -a')
