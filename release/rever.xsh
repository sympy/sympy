# -*-mode: python; flycheck-mode: nil -*-

$XONSH_SHOW_TRACEBACK = True
$RAISE_SUBPROC_ERROR = True

trace on

import os
import sys
import unicodedata
from collections import defaultdict, OrderedDict
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
from requests_oauthlib import OAuth2

from rever.activity import activity
from rever.conda import run_in_conda_env

cd ..

$ACTIVITIES = [
    # 'version_bump',
    '_version',
    'mailmap_update',
    'test_sympy',
    'source_tarball',
    'build_docs',
    'copy_release_files',
    'compare_tar_against_git',
    'test_tarball27',
    'test_tarball35',
    'test_tarball36',
    'test_tarball37',
    'print_authors',
    'md5',
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
    with run_in_conda_env(['python=3.6', 'mpmath']):
        ./bin/mailmap_update.py

@activity
def test_sympy():
    with run_in_conda_env(['mpmath', 'matplotlib>=2.2', 'numpy', 'scipy', 'theano',
        'ipython', 'gmpy2', 'fastcache', 'symengine', 'libgfortran', 'libgcc',
        'gcc', 'cython', 'tensorflow', 'llvmlite', 'wurlitzer', 'autowrap',
        'python-symengine=0.3.*', 'numexpr', 'antlr-python-runtime>=4.7,<4.8',
        'antlr>=4.7,<4.8'], 'sympy-tests'):

        ./setup.py test

@activity(deps={'_version', 'mailmap_update', 'test_sympy'})
def source_tarball():
    with run_in_conda_env(['mpmath', 'python=3.6'], 'sympy-release'):
        # Assumes this is run in Docker and git is already clean
        ./setup.py sdist --keep-temp

@activity(deps={'_version'})
def build_docs():
    with run_in_conda_env(['sphinx', 'docutils', 'numpy', 'mpmath'],
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


@activity(deps={'source_tarball', 'build_docs'})
def copy_release_files():
    ls dist
    cp dist/* /root/release/

@activity(deps={'source_tarball'})
def test_tarball27():
    test_tarball('2.7')

@activity(deps={'source_tarball'})
def test_tarball35():
    test_tarball('3.5')

@activity(deps={'source_tarball'})
def test_tarball36():
    test_tarball('3.6')

@activity(deps={'source_tarball'})
def test_tarball37():
    test_tarball('3.7')

@activity(deps={'source_tarball'})
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

@activity(deps={'source_tarball'})
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

@activity(deps={'source_tarball', 'build_docs'})
def md5():
    """
    Print the md5 sums of the release files
    """
    _md5(print_=True)

def _md5(print_=True, local=False):
    if local:
        out = $(md5sum @(release_files()))
    else:
        out = $(md5sum /root/release/*)
    # Remove the release/ part for printing. Useful for copy-pasting into the
    # release notes.
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out

@activity(deps={'mailmap_update', 'md5', 'print_authors', 'source_tarball', 'build_docs', 'compare_tar_against_git', 'test_tarball27', 'test_tarball35', 'test_tarball36', 'test_sympy'})
def release():
    pass

@activity(deps={'_version'})
def GitHub_release():
    _GitHub_release()

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
    if py_version not in {'2.7', '3.5', '3.6', '3.7'}: # TODO: Add win32
        raise ValueError("release must be one of 2.7, 3.5, 3.6, or 3.7 not %s" % py_version)


    with run_in_conda_env(['python=%s' % py_version], 'test-install-%s' % py_version):
        cp @('/root/release/{source}'.format(**tarball_format)) @("releasetar.tar".format(**tarball_format))
        tar xvf releasetar.tar

        cd @("/root/{source-orig-notar}".format(**tarball_format))
        python setup.py install
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

def release_files():
    """
    Returns the list of local release files
    """
    return glob.glob('release/release-' + $VERSION + '/*')

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
        ret = $(tar tf @("/root/release/{source}".format(**tarball_format)))
    elif file == 'html':
        ret = $(unzip -l @("/root/release/{html}".format(**tarball_format)))
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
    # Based on https://stackoverflow.com/a/13505966/161801
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

def _GitHub_release(username=None, user='sympy', token=None,
    token_file_path="~/.sympy/release-token", repo='sympy', draft=False):
    """
    Upload the release files to GitHub.

    The tag must be pushed up first. You can test on another repo by changing
    user and repo.
    """
    if not requests:
        error("requests and requests-oauthlib must be installed to upload to GitHub")

    release_text = GitHub_release_text()
    short_version = get_sympy_short_version()
    tag = 'sympy-' + version
    prerelease = short_version != version

    urls = URLs(user=user, repo=repo)
    if not username:
        username = input("GitHub username: ")
    token = load_token_file(token_file_path)
    if not token:
        username, password, token = GitHub_authenticate(urls, username, token)

    # If the tag in question is not pushed up yet, then GitHub will just
    # create it off of master automatically, which is not what we want.  We
    # could make it create it off the release branch, but even then, we would
    # not be sure that the correct commit is tagged.  So we require that the
    # tag exist first.
    if not check_tag_exists():
        sys.exit(red("The tag for this version has not been pushed yet. Cannot upload the release."))

    # See https://developer.github.com/v3/repos/releases/#create-a-release
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
        with open(os.path.join('release/release-' + $VERSION, tarball), 'rb') as f:
            result = query_GitHub(urls.release_uploads_url % release_id, username,
                password=None, token=token, data=f, params=params,
                headers=headers).json()

        print(green("Done"))

    # TODO: download the files and check that they have the right md5 sum

def _size(print_=True):
    """
    Print the sizes of the release files. Run locally.
    """
    out = $(du -h @(release_files()))
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out

def table():
    """
    Make an html table of the downloads.

    This is for pasting into the GitHub releases page. See GitHub_release().
    """
    tarball_formatter_dict = dict(tarball_format)
    shortversion = get_sympy_short_version()

    tarball_formatter_dict['version'] = shortversion

    md5s = [i.split('\t') for i in _md5(print_=False, local=True).split('\n')]
    md5s_dict = {name: md5 for md5, name in md5s}

    sizes = [i.split('\t') for i in _size(print_=False).split('\n')]
    sizes_dict = {name: size for size, name in sizes}

    table = []

    # https://docs.python.org/2/library/contextlib.html#contextlib.contextmanager. Not
    # recommended as a real way to generate html, but it works better than
    # anything else I've tried.
    @contextmanager
    def tag(name):
        table.append("<%s>" % name)
        yield
        table.append("</%s>" % name)
    @contextmanager
    def a_href(link):
        table.append("<a href=\"%s\">" % link)
        yield
        table.append("</a>")

    with tag('table'):
        with tag('tr'):
            for headname in ["Filename", "Description", "size", "md5"]:
                with tag("th"):
                    table.append(headname)

        for key in descriptions:
            name = get_tarball_name(key)
            with tag('tr'):
                with tag('td'):
                    with a_href('https://github.com/sympy/sympy/releases/download/sympy-%s/%s' % (version, name)):
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
files above.
"""
    out = out.format(shortversion=shortversion, htmltable=htmltable)
    print(blue("Here are the release notes to copy into the GitHub release "
        "Markdown form:"))
    print()
    print(out)
    return out

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
        username = input("Username: ")

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
        generate = input("> Generate API token? [Y/n] ")
        if generate.lower() in ["y", "ye", "yes", ""]:
            name = input("> Name of token on GitHub? [SymPy Release] ")
            if name == "":
                name = "SymPy Release"
            token = generate_token(urls, username, password, name=name)
            print("Your token is", token)
            print("Use this token from now on as GitHub_release:token=" + token +
                ",username=" + username)
            print(red("DO NOT share this token with anyone"))
            save = input("Do you want to save this token to a file [yes]? ")
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
    token_file = input("> Enter token file location [~/.sympy/release-token] ")
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
            OTP = input("Authentication code: ")
            return query_GitHub(url, username=username, password=password,
                token=token, data=data, OTP=OTP)

        raise AuthenticationFailed("invalid username or password")

    r.raise_for_status()
    return r

def check_tag_exists():
    """
    Check if the tag for this release has been uploaded yet.
    """
    tag = 'sympy-' + version
    all_tags = $(git ls-remote --tags origin)
    return tag in all_tags

descriptions = OrderedDict([
    ('source', "The SymPy source installer.",),
    ('html', '''Html documentation. This is the same as
the <a href="https://docs.sympy.org/latest/index.html">online documentation</a>.''',),
    ('pdf', '''Pdf version of the <a href="https://docs.sympy.org/latest/index.html"> html documentation</a>.''',),
    ])

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

    # See the README of the docs repo. We have to remove the old redirects,
    # move in the new docs, and create redirects.
    print("Removing redirects from previous version")
    rm -r @(previous_version)
    print("Moving previous latest docs to old version")
    mv latest @(previous_version)

    print("Unzipping docs into repo")
    unzip @(docs_zip) > /dev/null
    mv @(get_tarball_name('html-nozip')) @(version)

    print("Writing new version to releases.txt")
    with open(os.path.join(docs_location, "releases.txt"), 'a') as f:
        f.write("{version}:SymPy {version}\n".format(version=current_version))

    print("Generating indexes")
    ./generate_indexes.py
    mv @(current_version) latest

    print("Generating redirects")
    ./generate_redirects.py latest @(current_version)

    print("Committing")
    git add -A @(version) latest
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
    '.ci/blacklisted.json',
    '.editorconfig',
    # Code of conduct
    'CODE_OF_CONDUCT.md',
    # Pull request template
    'PULL_REQUEST_TEMPLATE.md',
    # Contributing guide
    'CONTRIBUTING.md',
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
    'bin/test_external_imports.py',
    'bin/test_executable.py',
    'bin/test_import',
    'bin/test_import.py',
    'bin/test_isolated',
    'bin/test_setup.py',
    'bin/test_travis.sh',
    # The notebooks are not ready for shipping yet. They need to be cleaned
    # up, and preferably doctested.  See also
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
    'examples/notebooks/Bezout_Dixon_resultant.ipynb',
    'examples/notebooks/IntegrationOverPolytopes.ipynb',
    'examples/notebooks/Macaulay_resultant.ipynb',
    'examples/notebooks/Sylvester_resultant.ipynb',
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
    'release/pull_and_run_rever.sh',
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
    'sympy.egg-info/entry_points.txt',
    }
