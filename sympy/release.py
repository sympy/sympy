from __future__ import annotations

__commit_id__: str | None
__version__: str

try:
    from sympy._version import __commit_id__, __version__
except ImportError:
    __commit_id__ = None
    __version__ = "0.0.dev0+uninstalled"
