#!/usr/bin/env python
"""Check the authors of PR commits against the repository's .mailmap.

Read one raw ``Name <email>`` identity per line from standard input. This
requires only the PR's commits, unlike the full-history mailmap_check.py.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


_IDENTITY = re.compile(r'[^<>]+ <[^<>]*>')
_CANONICAL = re.compile(r'^\s*([^<]+?)\s*<([^>]*)>')


def canonical_authors(mailmap):
    """Return the canonical identities recorded in .mailmap."""
    authors = set()
    for number, line in enumerate(mailmap.read_text(encoding='utf-8').splitlines(), 1):
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        match = _CANONICAL.match(line)
        if match is None:
            raise ValueError(f'Invalid .mailmap entry at line {number}')
        authors.add(f'{match[1].strip()} <{match[2]}>')
    return authors


def unrecognized_authors(authors, mailmap, expected):
    """Return (original, resolved) pairs not recognized by .mailmap."""
    authors = [author.strip() for author in authors]
    if len(authors) != expected:
        raise ValueError(
            f'Expected {expected} PR commits, received {len(authors)}. '
            'The commit list may be incomplete.')
    if any(not _IDENTITY.fullmatch(author) for author in authors):
        raise ValueError('PR contains an invalid commit author identity')

    result = subprocess.run(
        ['git', 'check-mailmap', '--stdin'],
        input='\n'.join(authors) + '\n', capture_output=True,
        text=True, encoding='utf-8', cwd=mailmap.parent, check=True)
    resolved = result.stdout.splitlines()
    if len(resolved) != len(authors):
        raise ValueError('git check-mailmap returned an unexpected number of authors')

    known = canonical_authors(mailmap)
    return [(original, canonical) for original, canonical in zip(authors, resolved)
            if canonical not in known]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected', required=True, type=int,
                        help='Number of commits reported by the PR event')
    parser.add_argument('--mailmap', type=Path,
                        default=Path(__file__).resolve().parent.parent / '.mailmap')
    args = parser.parse_args(argv)

    try:
        missing = unrecognized_authors(sys.stdin.read().splitlines(),
                                       args.mailmap, args.expected)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f'Author preflight failed: {exc}', file=sys.stderr)
        return 1

    if missing:
        print('PR commit authors not recognized by .mailmap:', file=sys.stderr)
        for original, canonical in missing:
            print(f'  {original} (resolved as {canonical})', file=sys.stderr)
        print('Add the author or missing alias to .mailmap.', file=sys.stderr)
        return 1

    print(f'All {args.expected} PR commit authors are recognized by .mailmap.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
