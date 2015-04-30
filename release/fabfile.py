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
*must* be changed in sympy/release.py for this to work).
"""
from __future__ import print_function

from collections import defaultdict, OrderedDict

from contextlib import contextmanager

from fabric.api import env, local, run, sudo, cd, hide, task
from fabric.contrib.files import exists
from fabric.colors import blue, red, green
from fabric.utils import error, warn

try:
    # Only works in newer versions of fabric
    env.colorize_errors = True
except AttributeError:
    pass

try:
    import requests
    from requests.auth import HTTPBasicAuth
    from requests_oauthlib import OAuth2
except ImportError:
    warn("requests and requests-oauthlib must be installed to upload to GitHub")
    requests = False

import unicodedata
import json
from getpass import getpass

import os
import stat
import sys

import time
import ConfigParser

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
    sudo("apt-get -y install graphviz inkscape texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra librsvg2-bin docbook2x")
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
    if os.path.exists("release"):
        error("release directory already exists locally. Remove it to continue.")

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
    Get the full version of SymPy being released (like 0.7.3.__weakrefoffset__rc1)
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
    parts = version.split('.')
    non_rc_parts = [i for i in parts if i.isdigit()]
    return '.'.join(non_rc_parts) # Remove any rc tags

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
    test_sympy()
    source_tarball()
    build_docs()
    copy_release_files()
    test_tarball('2')
    test_tarball('3')
    compare_tar_against_git()
    print_authors()

@task
def source_tarball():
    """
    Build the source tarball
    """
    with cd("/home/vagrant/repos/sympy"):
        run("git clean -dfx")
        run("./setup.py clean")
        run("./setup.py sdist --keep-temp")
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
                run("make man")
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
    # https://github.com/sympy/sympy/issues/7087.
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
    'bin/build_doc.sh',
    'bin/diagnose_imports',
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
    # This stuff :)
    'release/.gitignore',
    'release/README.md',
    'release/Vagrantfile',
    'release/fabfile.py',
    # This is just a distribute version of setup.py. Used mainly for setup.py
    # develop, which we don't care about in the release tarball
    'setupegg.py',
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
def size(file='*', print_=True):
    """
    Print the sizes of the release files
    """
    out = local("du -h release/" + file, capture=True)
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out

@task
def table():
    """
    Make an html table of the downloads.

    This is for pasting into the GitHub releases page. See GitHub_release().
    """
    # TODO: Add the file size
    tarball_formatter_dict = tarball_formatter()
    shortversion = get_sympy_short_version()

    tarball_formatter_dict['version'] = shortversion

    md5s = [i.split('\t') for i in md5(print_=False).split('\n')]
    md5s_dict = {name: md5 for md5, name in md5s}

    sizes = [i.split('\t') for i in size(print_=False).split('\n')]
    sizes_dict = {name: size for size, name in sizes}

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
            for headname in ["Filename", "Description", "size", "md5"]:
                with tag("th"):
                    table.append(headname)

        for key in descriptions:
            name = get_tarball_name(key)
            with tag('tr'):
                with tag('td'):
                    with tag('b'):
                        table.append(name)
                with tag('td'):
                    table.append(descriptions[key].format(**tarball_formatter_dict))
                with tag('td'):
                    table.append(sizes_dict[name])
                with tag('td'):
                    table.append(md5s_dict[name])

    out = ' '.join(table)
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
                # If the tagged commit *is* a merge commit, just comment this
                # out, and make sure `fab vagrant get_previous_version_tag` is correct
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

@task
def check_tag_exists():
    """
    Check if the tag for this release has been uploaded yet.
    """
    version = get_sympy_version()
    tag = 'sympy-' + version
    with cd("/home/vagrant/repos/sympy"):
        all_tags = run("git ls-remote --tags origin")
    return tag in all_tags

# ------------------------------------------------
# Updating websites

@task
def update_websites():
    """
    Update various websites owned by SymPy.

    So far, supports the docs and sympy.org
    """
    update_docs()
    update_sympy_org()

def get_location(location):
    """
    Read/save a location from the configuration file.
    """
    locations_file = os.path.expanduser('~/.sympy/sympy-locations')
    config = ConfigParser.SafeConfigParser()
    config.read(locations_file)
    the_location = config.has_option("Locations", location) and config.get("Locations", location)
    if not the_location:
        the_location = raw_input("Where is the SymPy {location} directory? ".format(location=location))
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

@task
def update_docs(docs_location=None):
    """
    Update the docs hosted at docs.sympy.org
    """
    docs_location = docs_location or get_location("docs")

    print("Docs location:", docs_location)

    # Check that the docs directory is clean
    local("cd {docs_location} && git diff --exit-code > /dev/null".format(docs_location=docs_location))
    local("cd {docs_location} && git diff --cached --exit-code > /dev/null".format(docs_location=docs_location))

    # See the README of the docs repo. We have to remove the old redirects,
    # move in the new docs, and create redirects.
    current_version = get_sympy_version()
    previous_version = get_previous_version_tag().lstrip('sympy-')
    print("Removing redirects from previous version")
    local("cd {docs_location} && rm -r {previous_version}".format(docs_location=docs_location,
        previous_version=previous_version))
    print("Moving previous latest docs to old version")
    local("cd {docs_location} && mv latest {previous_version}".format(docs_location=docs_location,
        previous_version=previous_version))

    print("Unzipping docs into repo")
    release_dir = os.path.abspath(os.path.expanduser(os.path.join(os.path.curdir, 'release')))
    docs_zip = os.path.abspath(os.path.join(release_dir, get_tarball_name('html')))
    local("cd {docs_location} && unzip {docs_zip} > /dev/null".format(docs_location=docs_location,
        docs_zip=docs_zip))
    local("cd {docs_location} && mv {docs_zip_name} {version}".format(docs_location=docs_location,
        docs_zip_name=get_tarball_name("html-nozip"), version=current_version))

    print("Writing new version to releases.txt")
    with open(os.path.join(docs_location, "releases.txt"), 'a') as f:
        f.write("{version}:SymPy {version}\n".format(version=current_version))

    print("Generating indexes")
    local("cd {docs_location} && ./generate_indexes.py".format(docs_location=docs_location))
    local("cd {docs_location} && mv {version} latest".format(docs_location=docs_location,
        version=current_version))

    print("Generating redirects")
    local("cd {docs_location} && ./generate_redirects.py latest {version} ".format(docs_location=docs_location,
        version=current_version))

    print("Committing")
    local("cd {docs_location} && git add -A {version} latest".format(docs_location=docs_location,
        version=current_version))
    local("cd {docs_location} && git commit -a -m \'Updating docs to {version}\'".format(docs_location=docs_location,
        version=current_version))

    print("Pushing")
    local("cd {docs_location} && git push origin".format(docs_location=docs_location))

@task
def update_sympy_org(website_location=None):
    """
    Update sympy.org

    This just means adding an entry to the news section.
    """
    website_location = website_location or get_location("sympy.github.com")

    # Check that the website directory is clean
    local("cd {website_location} && git diff --exit-code > /dev/null".format(website_location=website_location))
    local("cd {website_location} && git diff --cached --exit-code > /dev/null".format(website_location=website_location))

    release_date = time.gmtime(os.path.getctime(os.path.join("release",
        tarball_formatter()['source'])))
    release_year = str(release_date.tm_year)
    release_month = str(release_date.tm_mon)
    release_day = str(release_date.tm_mday)
    version = get_sympy_version()

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
    local("cd {website_location} && ./generate".format(website_location=website_location))

    print("Committing")
    local("cd {website_location} && git commit -a -m \'Add {version} to the news\'".format(website_location=website_location,
        version=version))

    print("Pushing")
    local("cd {website_location} && git push origin".format(website_location=website_location))

# ------------------------------------------------
# Uploading

@task
def upload():
    """
    Upload the files everywhere (PyPI and GitHub)

    """
    distutils_check()
    GitHub_release()
    pypi_register()
    pypi_upload()
    test_pypi(2)
    test_pypi(3)

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
    Upload files to PyPI. You will need to enter a password.
    """
    with cd("/home/vagrant/repos/sympy"):
        # See http://stackoverflow.com/a/17657183/161801
        run("python setup.py sdist --dry-run upload")

@task
def test_pypi(release='2'):
    """
    Test that the sympy can be pip installed, and that sympy imports in the
    install.
    """
    # This function is similar to test_tarball()

    version = get_sympy_version()

    release = str(release)

    if release not in {'2', '3'}: # TODO: Add win32
        raise ValueError("release must be one of '2', '3', not %s" % release)

    venv = "/home/vagrant/repos/test-{release}-pip-virtualenv".format(release=release)

    with use_venv(release):
        make_virtualenv(venv)
        with virtualenv(venv):
            run("pip install sympy")
            run('python -c "import sympy; assert sympy.__version__ == \'{version}\'"'.format(version=version))

@task
def GitHub_release_text():
    """
    Generate text to put in the GitHub release Markdown box
    """
    shortversion = get_sympy_short_version()
    htmltable = table()
    out = """\
See https://github.com/sympy/sympy/wiki/release-notes-for-{shortversion} for the release notes.

{htmltable}

**Note**: Do not download the **Source code (zip)** or the **Source code (tar.gz)**
files below.
"""
    out = out.format(shortversion=shortversion, htmltable=htmltable)
    print(blue("Here are the release notes to copy into the GitHub release "
        "Markdown form:", bold=True))
    print()
    print(out)
    return out

@task
def GitHub_release(username=None, user='sympy', token=None,
    token_file_path="~/.sympy/release-token", repo='sympy', draft=False):
    """
    Upload the release files to GitHub.

    The tag must be pushed up first. You can test on another repo by changing
    user and repo.
    """
    if not requests:
        error("requests and requests-oauthlib must be installed to upload to GitHub")

    release_text = GitHub_release_text()
    version = get_sympy_version()
    short_version = get_sympy_short_version()
    tag = 'sympy-' + version
    prerelease = short_version != version

    urls = URLs(user=user, repo=repo)
    if not username:
        username = raw_input("GitHub username: ")
    token = load_token_file(token_file_path)
    if not token:
        username, password, token = GitHub_authenticate(urls, username, token)

    # If the tag in question is not pushed up yet, then GitHub will just
    # create it off of master automatically, which is not what we want.  We
    # could make it create it off the release branch, but even then, we would
    # not be sure that the correct commit is tagged.  So we require that the
    # tag exist first.
    if not check_tag_exists():
        error("The tag for this version has not been pushed yet. Cannot upload the release.")

    # See http://developer.github.com/v3/repos/releases/#create-a-release
    # First, create the release
    post = {}
    post['tag_name'] = tag
    post['name'] = "SymPy " + version
    post['body'] = release_text
    post['draft'] = draft
    post['prerelease'] = prerelease

    print("Creating release for tag", tag, end=' ')

    result = query_GitHub(urls.releases_url, username, password=None,
        token=token, data=json.dumps(post)).json()
    release_id = result['id']

    print(green("Done"))

    # Then, upload all the files to it.
    for key in descriptions:
        tarball = get_tarball_name(key)

        params = {}
        params['name'] = tarball

        if tarball.endswith('gz'):
            headers = {'Content-Type':'application/gzip'}
        elif tarball.endswith('pdf'):
            headers = {'Content-Type':'application/pdf'}
        elif tarball.endswith('zip'):
            headers = {'Content-Type':'application/zip'}
        else:
            headers = {'Content-Type':'application/octet-stream'}

        print("Uploading", tarball, end=' ')
        sys.stdout.flush()
        with open(os.path.join("release", tarball), 'rb') as f:
            result = query_GitHub(urls.release_uploads_url % release_id, username,
                password=None, token=token, data=f, params=params,
                headers=headers).json()

        print(green("Done"))

    # TODO: download the files and check that they have the right md5 sum

def GitHub_check_authentication(urls, username, password, token):
    """
    Checks that username & password is valid.
    """
    query_GitHub(urls.api_url, username, password, token)

def GitHub_authenticate(urls, username, token=None):
    _login_message = """\
Enter your GitHub username & password or press ^C to quit. The password
will be kept as a Python variable as long as this script is running and
https to authenticate with GitHub, otherwise not saved anywhere else:\
"""
    if username:
        print("> Authenticating as %s" % username)
    else:
        print(_login_message)
        username = raw_input("Username: ")

    authenticated = False

    if token:
        print("> Authenticating using token")
        try:
            GitHub_check_authentication(urls, username, None, token)
        except AuthenticationFailed:
            print(">     Authentication failed")
        else:
            print(">     OK")
            password = None
            authenticated = True

    while not authenticated:
        password = getpass("Password: ")
        try:
            print("> Checking username and password ...")
            GitHub_check_authentication(urls, username, password, None)
        except AuthenticationFailed:
            print(">     Authentication failed")
        else:
            print(">     OK.")
            authenticated = True

    if password:
        generate = raw_input("> Generate API token? [Y/n] ")
        if generate.lower() in ["y", "ye", "yes", ""]:
            name = raw_input("> Name of token on GitHub? [SymPy Release] ")
            if name == "":
                name = "SymPy Release"
            token = generate_token(urls, username, password, name=name)
            print("Your token is", token)
            print("Use this token from now on as GitHub_release:token=" + token +
                ",username=" + username)
            print(red("DO NOT share this token with anyone"))
            save = raw_input("Do you want to save this token to a file [yes]? ")
            if save.lower().strip() in ['y', 'yes', 'ye', '']:
                save_token_file(token)

    return username, password, token

def generate_token(urls, username, password, OTP=None, name="SymPy Release"):
    enc_data = json.dumps(
        {
            "scopes": ["public_repo"],
            "note": name
        }
    )

    url = urls.authorize_url
    rep = query_GitHub(url, username=username, password=password,
        data=enc_data).json()
    return rep["token"]

def save_token_file(token):
    token_file = raw_input("> Enter token file location [~/.sympy/release-token] ")
    token_file = token_file or "~/.sympy/release-token"

    token_file_expand = os.path.expanduser(token_file)
    token_file_expand = os.path.abspath(token_file_expand)
    token_folder, _ = os.path.split(token_file_expand)

    try:
        if not os.path.isdir(token_folder):
            os.mkdir(token_folder, 0o700)
        with open(token_file_expand, 'w') as f:
            f.write(token + '\n')
        os.chmod(token_file_expand, stat.S_IREAD | stat.S_IWRITE)
    except OSError as e:
        print("> Unable to create folder for token file: ", e)
        return
    except IOError as e:
        print("> Unable to save token file: ", e)
        return

    return token_file

def load_token_file(path="~/.sympy/release-token"):
    print("> Using token file %s" % path)

    path = os.path.expanduser(path)
    path = os.path.abspath(path)

    if os.path.isfile(path):
        try:
            with open(path) as f:
                token = f.readline()
        except IOError:
            print("> Unable to read token file")
            return
    else:
        print("> Token file does not exist")
        return

    return token.strip()

class URLs(object):
    """
    This class contains URLs and templates which used in requests to GitHub API
    """

    def __init__(self, user="sympy", repo="sympy",
        api_url="https://api.github.com",
        authorize_url="https://api.github.com/authorizations",
        uploads_url='https://uploads.github.com',
        main_url='https://github.com'):
        """Generates all URLs and templates"""

        self.user = user
        self.repo = repo
        self.api_url = api_url
        self.authorize_url = authorize_url
        self.uploads_url = uploads_url
        self.main_url = main_url

        self.pull_list_url = api_url + "/repos" + "/" + user + "/" + repo + "/pulls"
        self.issue_list_url = api_url + "/repos/" + user + "/" + repo + "/issues"
        self.releases_url = api_url + "/repos/" + user + "/" + repo + "/releases"
        self.single_issue_template = self.issue_list_url + "/%d"
        self.single_pull_template = self.pull_list_url + "/%d"
        self.user_info_template = api_url + "/users/%s"
        self.user_repos_template = api_url + "/users/%s/repos"
        self.issue_comment_template = (api_url + "/repos" + "/" + user + "/" + repo + "/issues/%d" +
            "/comments")
        self.release_uploads_url = (uploads_url + "/repos/" + user + "/" +
            repo + "/releases/%d" + "/assets")
        self.release_download_url = (main_url + "/" + user + "/" + repo +
            "/releases/download/%s/%s")


class AuthenticationFailed(Exception):
    pass

def query_GitHub(url, username=None, password=None, token=None, data=None,
    OTP=None, headers=None, params=None, files=None):
    """
    Query GitHub API.

    In case of a multipage result, DOES NOT query the next page.

    """
    headers = headers or {}

    if OTP:
        headers['X-GitHub-OTP'] = OTP

    if token:
        auth = OAuth2(client_id=username, token=dict(access_token=token,
            token_type='bearer'))
    else:
        auth = HTTPBasicAuth(username, password)
    if data:
        r = requests.post(url, auth=auth, data=data, headers=headers,
            params=params, files=files)
    else:
        r = requests.get(url, auth=auth, headers=headers, params=params, stream=True)

    if r.status_code == 401:
        two_factor = r.headers.get('X-GitHub-OTP')
        if two_factor:
            print("A two-factor authentication code is required:", two_factor.split(';')[1].strip())
            OTP = raw_input("Authentication code: ")
            return query_GitHub(url, username=username, password=password,
                token=token, data=data, OTP=OTP)

        raise AuthenticationFailed("invalid username or password")

    r.raise_for_status()
    return r

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
