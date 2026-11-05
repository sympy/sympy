"""Regression tests for the lightweight PR author preflight."""

import subprocess
import sys
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[3] / 'bin' / 'check_pr_authors.py'


def test_check_pr_authors(tmp_path):
    subprocess.run(['git', 'init', '-q'], cwd=tmp_path, check=True)
    (tmp_path / '.mailmap').write_text(
        'Alice Smith <alice@example.com>\n'
        'Alice Smith <alice@example.com> Alice <old@example.com>\n'
        'Bob Jones <bob@example.com>\n', encoding='utf-8')

    def check(authors, expected):
        return subprocess.run(
            [sys.executable, str(_SCRIPT), '--expected', str(expected),
             '--mailmap', str(tmp_path / '.mailmap')],
            input='\n'.join(authors) + '\n', capture_output=True,
            text=True, encoding='utf-8', cwd=tmp_path)

    # Existing canonical identities, including identities with a known alias.
    ok = check(['Alice Smith <alice@example.com>',
                'Alice <old@example.com>'], 2)
    assert ok.returncode == 0, ok.stderr

    # A new contributor is accepted as soon as their entry is added.
    ok = check(['Bob Jones <bob@example.com>'], 1)
    assert ok.returncode == 0, ok.stderr

    # An existing author using an unlisted email still needs a new alias.
    bad = check(['Alice Smith <new@example.com>'], 1)
    assert bad.returncode == 1
    assert 'new@example.com' in bad.stderr

    # Don't silently succeed if the API returned an incomplete commit list.
    bad = check(['Alice Smith <alice@example.com>'], 2)
    assert bad.returncode == 1
    assert 'incomplete' in bad.stderr

    # Empty results must not accidentally pass an expected nonempty PR.
    assert check([], 1).returncode == 1
