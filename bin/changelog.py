import os
import sys
import re
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from git import Repo, Actor
from collections import defaultdict

HEADERS = {
    'M': 'Major changes',
    'b': 'Backwards compatibility breaks and deprecations',
    'n': 'New features',
    'm': 'Minor changes'
}
PREFIX = '* '
SUFFIX = ' (#%s)\n' # % pull request number

def get_build_information():
    """
    Return (whether it is a push build, PR number)
    """
    event_type = os.environ['TRAVIS_EVENT_TYPE']
    branch = os.environ['TRAVIS_BRANCH']
    if branch != 'master':
        print('Target branch is not master, skipping changelog')
        sys.exit()
    elif event_type == 'push':
        match = re.match(r'Merge pull request #(\d+)',
            os.environ['TRAVIS_COMMIT_MESSAGE']);
        if not match:
            sys.exit('Cannot find pull request number!')
        return (True, match.group(1))
    elif event_type == 'pull_request':
        return (False, os.environ['TRAVIS_PULL_REQUEST'])
    elif event_type == 'cron' or event_type == 'api':
        print('Not a push or PR build, skipping changelog')
        sys.exit()
    else:
        sys.exit('Unknown event type!')

def get_pr_desc(pr_number):
    """
    Return pull request description retrieved using GitHub API
    """
    s = requests.Session()
    retry = Retry(total=5, read=5, connect=5, backoff_factor=0.1,
        status_forcelist=(500, 502, 503, 504))
    adapter = HTTPAdapter(max_retries=retry)
    s.mount('https://', adapter)
    r = s.get('https://api.github.com/repos/' +
        os.environ['TRAVIS_REPO_SLUG'] + '/pulls/' + pr_number)
    r.raise_for_status()
    return r.json()['body']

def get_changelog(data):
    """
    Return changelogs parsed from string
    """
    changelogs = defaultdict(list)
    start = False
    for line in data.splitlines():
        if 'This PR changes' in line:
            start = True
        if start:
            if ':' in line:
                type, message = line.split(':', 1)
                changelogs[type.strip()].append(message.strip())
            elif '[CHANGELOG END]' in line:
                break
            elif '[skip changelog]' in line:
                print('[skip changelog] found, skipping changelog')
                sys.exit()
    else:
        sys.exit('Changelog not detected! Please use pull request template.')
    if all(len(changelogs[t]) == 0 for t in HEADERS):
        sys.exit('Changelog not found! Please add a changelog.')
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
        if line.startswith('##'):
            # last heading is assumed to be not a changelog header
            if current:
                suffix = SUFFIX % pr_number
                entry = (PREFIX + (suffix + PREFIX).join(changelogs[header]) +
                    suffix + '\n')
                if not is_prev_empty:
                    contents[i] = contents[i] + entry
                else:
                    contents[n] = entry
            for header in HEADERS:
                if HEADERS[header] in line:
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
    update, pr_number = get_build_information()
    pr_desc = get_pr_desc(pr_number)
    changelogs = get_changelog(pr_desc)

    if True:
        # release-notes-bot:' + os.environ['BOT_TOKEN'] + '@
        # TODO : encrypt token
        repo = Repo.clone_from('https://github.com/' + os.environ['TRAVIS_REPO_SLUG'] +
            '.wiki.git', os.path.abspath('wiki'))

        rel_notes_path = os.path.join(repo.working_tree_dir,
            get_release_notes_filename())
        update_release_notes(rel_notes_path, changelogs, pr_number)

        with open(rel_notes_path, 'r') as f:
            print(f.read(), end="")
        #repo.index.add([rel_notes_path])
        #bot = Actor('Release Notes Bot', '33476835+release-notes-bot@users.noreply.github.com')
        #repo.index.commit('Update release notes after pull request #%s' %
        #    pr_number, author=bot, committer=bot)
        #repo.remotes.origin.push()
