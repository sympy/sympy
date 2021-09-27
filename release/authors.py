#!/usr/bin/env python3

import os
from pathlib import Path
from subprocess import check_output
import unicodedata


def main(version, prevversion, outdir):
    """
    Print authors text to put at the bottom of the release notes
    """
    outdir = Path(outdir)
    authors, authorcount, newauthorcount = get_authors(version, prevversion)

    authors_text = f"""## Authors

The following people contributed at least one patch to this release (names are
given in alphabetical order by last name). A total of {authorcount} people
contributed to this release. People with a * by their names contributed a
patch for the first time for this release; {newauthorcount} people contributed
for the first time for this release.

Thanks to everyone who contributed to this release!
"""

    authors_lines = []
    for name in authors:
        authors_lines.append("- " + name)

    authors_text += '\n'.join(authors_lines)

    # Output to file and to screen
    with open(outdir / 'authors.txt', 'w') as authorsfile:
        authorsfile.write(authors_text)

    print()
    print(blue("Here are the authors to put at the bottom of the release notes."))
    print()
    print(authors_text)


def blue(text):
    return "\033[34m%s\033[0m" % text


def get_authors(version, prevversion):
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

    # The get_previous_version function can be flakey so we require the
    # previous version to be provided explicitly by the caller.
    #
    #old_release_tag = get_previous_version_tag(version)
    old_release_tag = 'sympy-' + prevversion

    out = check_output(['git', '--no-pager', 'log', old_release_tag + '..', '--format=%aN'])
    releaseauthors = set(out.decode('utf-8').strip().split('\n'))
    out = check_output(['git', '--no-pager', 'log', old_release_tag, '--format=%aN'])
    priorauthors = set(out.decode('utf-8').strip().split('\n'))

    releaseauthors = {name.strip() for name in releaseauthors if name.strip()}
    priorauthors = {name.strip() for name in priorauthors if name.strip()}
    newauthors = releaseauthors - priorauthors
    starred_newauthors = {name + "*" for name in newauthors}
    authors = releaseauthors - newauthors | starred_newauthors
    return (sorted(authors, key=lastnamekey), len(releaseauthors), len(newauthors))


def get_previous_version_tag(version):
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
    shortversion = get_sympy_short_version(version)
    curcommit = "HEAD"
    while True:
        cmdline = f'git describe --abbrev=0 --tags {curcommit}'
        print(cmdline)
        curtag = check_output(cmdline.split()).decode('utf-8').strip()
        if shortversion in curtag:
            # If the tagged commit is a merge commit, we cannot be sure
            # that it will go back in the right direction. This almost
            # never happens, so just error
            cmdline = f'git rev-list --parents -n 1 {curtag}'
            print(cmdline)
            parents = check_output(cmdline.split()).decode('utf-8').strip().split()
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


def get_sympy_short_version(version):
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


if __name__ == "__main__":
    import sys
    sys.exit(main(*sys.argv[1:]))
