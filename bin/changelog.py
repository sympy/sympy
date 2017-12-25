#!/usr/bin/env python
import os
import sys
import re
import requests
import subprocess
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from collections import defaultdict

HEADERS = [
    'Major changes',
    'Backwards compatibility breaks',
    'Minor changes'
]
PREFIX = '* '
SUFFIX = ' ([#{pr_number}](../pull/{pr_number}))\n'

def red(text):
    return "\033[31m%s\033[0m" % text

def green(text):
    return "\033[32m%s\033[0m" % text

def blue(text):
    return "\033[34m%s\033[0m" % text

def get_build_information():
    """
    Check event type and return PR number
    """
    event_type = os.environ['TRAVIS_EVENT_TYPE']
    if event_type == 'push':
        match = re.match(r'Merge pull request #(\d+)',
            os.environ['TRAVIS_COMMIT_MESSAGE']);
        if not match:
            sys.exit(red('Cannot find pull request number!'))
        return match.group(1)
    elif event_type == 'pull_request':
        return os.environ['TRAVIS_PULL_REQUEST']
    elif event_type == 'cron' or event_type == 'api':
        print(green('Not a push or PR build, skipping changelog'))
        sys.exit()
    else:
        sys.exit(red('Unknown event type!'))

def request_https_get(url):
    """
    Make HTTPS GET request and return response
    """
    s = requests.Session()
    retry = Retry(total=5, read=5, connect=5, backoff_factor=0.1,
        status_forcelist=(500, 502, 503, 504))
    adapter = HTTPAdapter(max_retries=retry)
    s.mount('https://', adapter)
    r = s.get(url)
    r.raise_for_status()
    return r

def get_pr_desc(pr_number):
    """
    Retrieve pull request description using GitHub API
    """
    r = request_https_get('https://api.github.com/repos/' +
        os.environ['TRAVIS_REPO_SLUG'] + '/pulls/' + pr_number)
    return r.json()['body']

def get_changelog(data):
    """
    Parse changelogs from a string
    """
    changelogs = defaultdict(list)
    current = None
    for line in data.splitlines():
        if current:
            if 'Add entry(ies)' in line:
                _, answer = line.split('?', 1)
                answer = re.sub('[^a-zA-Z]+', '', answer).lower()
                if answer in ['no', 'n', 'false', 'f']:
                    print(green('Skipping changelog'))
                    sys.exit()
                break
            elif line.startswith('*'):
                _, message = line.split('*', 1)
                message = message.strip()
                if message:
                    changelogs[header].append(message)
            elif line.startswith('#### '):
                for header in HEADERS:
                    if header in line:
                        current = header
                        break
        elif 'Brief description' in line:
            current = '_'
    else:
        sys.exit(red('Changelog not detected! Please use pull request template.'))
    count = sum(len(changelogs[t]) for t in HEADERS)
    if count == 0:
        sys.exit(red('Changelog not found! Please add a changelog.'))
    print(green('%s changelog entr%s found!' % (count, 'y' if count == 1 else 'ies')))
    return changelogs

def get_release_notes_filename():
    """
    Return filename of release notes for current development version
    """
    import sympy
    v = re.match(r'\d+(?:(?:\.\d+)*(?:\.[1-9]\d*)|\.0)', sympy.__version__).group()
    return 'Release-Notes-for-' + v + '.md'

def update_release_notes(rel_notes_path, changelogs, pr_number):
    """
    Update release notes
    """
    rel_notes = open(rel_notes_path, 'r+')
    contents = rel_notes.readlines()

    current = None
    for i in range(len(contents)):
        line = contents[i].strip()
        is_empty = (not line or line == '*')
        if line.startswith('## '):
            # last heading is assumed to be not a changelog header
            if changelogs[current]:
                suffix = SUFFIX.format(pr_number=pr_number)
                entry = (PREFIX + (suffix + PREFIX).join(changelogs[current]) +
                    suffix + '\n')
                if not is_prev_empty:
                    contents[i] = contents[i] + entry
                else:
                    contents[n] = entry
            for header in HEADERS:
                if header in line:
                    current = header
                    break
        elif current and not is_prev_empty and is_empty:
            n = i
        is_prev_empty = is_empty

    rel_notes.seek(0)
    rel_notes.writelines(contents)
    rel_notes.truncate()
    rel_notes.close()

if __name__ == '__main__':
    ON_TRAVIS = os.environ.get('TRAVIS', 'false') == 'true'

    if ON_TRAVIS:
        pr_number = get_build_information()
    else:
        while True:
            test_repo = input("Enter GitHub repository name: ")
            if test_repo.count("/") == 1:
                os.environ['TRAVIS_REPO_SLUG'] = test_repo
                break
        while True:
            pr_number = input("Enter PR number: ")
            if pr_number.isdigit():
                break

    pr_desc = get_pr_desc(pr_number)
    changelogs = get_changelog(pr_desc)

    rel_notes_name = get_release_notes_filename()
    rel_notes_path = os.path.abspath(rel_notes_name)

    if not ON_TRAVIS:
        print(green('Parsed changelogs:'))
        print(str(changelogs))
        print(green('Downloading %s' % rel_notes_name))
        r = request_https_get('https://raw.githubusercontent.com/wiki/' +
            os.environ['TRAVIS_REPO_SLUG'] + '/' + rel_notes_name)
        with open(rel_notes_path, 'w') as f:
            f.write(r.text)

    print(green('Updating %s' % rel_notes_path))
    update_release_notes(rel_notes_path, changelogs, pr_number)

    if ON_TRAVIS:
        print(green('Staging updated release notes'))
        command = ['git', 'add', rel_notes_path]
        print(blue(' '.join(command)))
        subprocess.run(command, check=True)
