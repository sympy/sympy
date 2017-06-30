# -*-mode: python; flycheck-mode: nil -*-

$XONSH_SHOW_TRACEBACK = True

from rever.activity import activity

cd ..

$ACTIVITIES = [
    # 'version_bump',
    'mailmap_update',
    # 'tag',
]

$TAG_PUSH = False

$VERSION_BUMP_PATTERNS = [
    ('sympy/release.py', r'__version__ = ".*"', r'__version__ = "$VERSION"'),
    ]

@activity
def mailmap_update():
    ./bin/mailmap_update.py
