# -*-mode: python; flycheck-mode: nil -*-

$ACTIVITIES = {
    'version_bump',
    'tag',
}

$TAG_PUSH = False

$VERSION_BUMP_PATTERNS = [
    ('sympy/release.py', r'__version__ = ".*"', r'__version__ = "$VERSION"'),
    ]
