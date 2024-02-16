#!/usr/bin/env python

import os
import json
from subprocess import check_output
from collections import OrderedDict, defaultdict
from collections.abc import Mapping
import glob
from contextlib import contextmanager

import requests
from requests_oauthlib import OAuth2


def main(version, push=None):
    """
    WARNING: If push is given as --push then this will push the release to
    github.
    """
    push = push == '--push'
    _GitHub_release(version, push)


def error(msg):
    raise ValueError(msg)


def blue(text):
    return "\033[34m%s\033[0m" % text

def red(text):
    return "\033[31m%s\033[0m" % text

def green(text):
    return "\033[32m%s\033[0m" % text


def _GitHub_release(version, push, username=None, user='sympy', token=None,
    token_file_path="~/.sympy/release-token", repo='sympy', draft=False):
    """
    Upload the release files to GitHub.

    The tag must be pushed up first. You can test on another repo by changing
    user and repo.
    """
    if not requests:
        error("requests and requests-oauthlib must be installed to upload to GitHub")

    release_text = GitHub_release_text(version)
    short_version = get_sympy_short_version(version)
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
    if not check_tag_exists(version):
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

    if push:
        result = query_GitHub(urls.releases_url, username, password=None,
            token=token, data=json.dumps(post)).json()
        release_id = result['id']
    else:
        print(green("Not pushing!"))

    print(green("Done"))

    # Then, upload all the files to it.
    for key in descriptions:
        tarball = get_tarball_name(key, version)

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
        with open(os.path.join('release/release-' + version, tarball), 'rb') as f:
            if push:
                result = query_GitHub(urls.release_uploads_url % release_id, username,
                    password=None, token=token, data=f, params=params,
                    headers=headers).json()
            else:
                print(green("Not uploading!"))

        print(green("Done"))

    # TODO: download the files and check that they have the right sha256 sum


def GitHub_release_text(version):
    """
    Generate text to put in the GitHub release Markdown box
    """
    shortversion = get_sympy_short_version(version)
    htmltable = table(version)
    out = """\
See https://github.com/sympy/sympy/wiki/release-notes-for-{shortversion} for the release notes.

{htmltable}

**Note**: Do not download the **Source code (zip)** or the **Source code (tar.gz)**
files below.
"""
    out = out.format(shortversion=shortversion, htmltable=htmltable)
    print(blue("Here are the release notes to copy into the GitHub release "
        "Markdown form:"))
    print()
    print(out)
    return out


def get_sympy_short_version(version):
    """
    Get the short version of SymPy being released, not including any rc tags
    (like 0.7.3)
    """
    parts = version.split('.')
    # Remove rc tags e.g. 1.10rc1 -> [1, 10]
    lastpart = ''
    for dig in parts[-1]:
        if dig.isdigit():
            lastpart += dig
        else:
            break
    parts[-1] = lastpart
    return '.'.join(parts)


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


def run(*cmdline, cwd=None):
    """
    Run command in subprocess and get lines of output
    """
    return check_output(cmdline, encoding='utf-8', cwd=cwd).splitlines()


def check_tag_exists(version):
    """
    Check if the tag for this release has been uploaded yet.
    """
    tag = 'sympy-' + version
    all_tag_lines = run('git', 'ls-remote', '--tags', 'origin')
    return any(tag in tag_line for tag_line in all_tag_lines)


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


def GitHub_check_authentication(urls, username, password, token):
    """
    Checks that username & password is valid.
    """
    query_GitHub(urls.api_url, username, password, token)


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
        auth = OAuth2(client_id=username, token={"access_token": token,
            "token_type": 'bearer'})
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


def table(version):
    """
    Make an html table of the downloads.

    This is for pasting into the GitHub releases page. See GitHub_release().
    """
    tarball_formatter_dict = dict(_tarball_format(version))
    shortversion = get_sympy_short_version(version)

    tarball_formatter_dict['version'] = shortversion

    sha256s = [i.split('\t') for i in _sha256(version, print_=False, local=True).split('\n')]
    sha256s_dict = {name: sha256 for sha256, name in sha256s}

    sizes = [i.split('\t') for i in _size(version, print_=False).split('\n')]
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
            for headname in ["Filename", "Description", "size", "sha256"]:
                with tag("th"):
                    table.append(headname)

        for key in descriptions:
            name = get_tarball_name(key, version)
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
                    table.append(sha256s_dict[name])

    out = ' '.join(table)
    return out

descriptions = OrderedDict([
    ('source', "The SymPy source installer.",),
    ('wheel', "A wheel of the package.",),
    ('html', '''Html documentation. This is the same as
the <a href="https://docs.sympy.org/latest/index.html">online documentation</a>.''',),
    ('pdf', '''Pdf version of the <a href="https://docs.sympy.org/latest/index.html"> html documentation</a>.''',),
    ])


def _size(version, print_=True):
    """
    Print the sizes of the release files. Run locally.
    """
    out = run(*(['du', '-h'] + release_files(version)))
    out = [i.split() for i in out]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out


def _sha256(version, print_=True, local=False):
    if local:
        out = run(*(['shasum', '-a', '256'] + release_files(version)))
    else:
        raise ValueError('Should not get here...')
        # out = run(*(['shasum', '-a', '256', '/root/release/*']))
    # Remove the release/ part for printing. Useful for copy-pasting into the
    # release notes.
    out = [i.split() for i in out]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])
    if print_:
        print(out)
    return out


def get_tarball_name(file, version):
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


def release_files(version):
    """
    Returns the list of local release files
    """
    paths = glob.glob('release/release-' + version + '/*')
    if not paths:
        raise ValueError("No release files found")
    return paths


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

    def __init__(self, version):
        self.version = version

    def __getitem__(self, name):
        return get_tarball_name(name, self.version)

    def __iter__(self):
        return iter(tarball_name_types)

    def __len__(self):
        return len(tarball_name_types)


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
