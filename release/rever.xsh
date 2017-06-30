# -*-mode: python; flycheck-mode: nil -*-

from rever.activity import activity

$ACTIVITIES = [
    'version_bump',
    'mailmap_update',
    # 'tag',
]

$TAG_PUSH = False

$VERSION_BUMP_PATTERNS = [
    ('sympy/release.py', r'__version__ = ".*"', r'__version__ = "$VERSION"'),
    ]

@activity
def mailmap_update():
    ../bin/mailmap_update.py
