def get_intersphinx_mapping(packages=None):
    """Get intersphinx mapping for specified packages."""
    mapping = {
        "matplotlib": ("https://matplotlib.org/stable/", None),
        "mpmath": ("https://mpmath.org/doc/current/", None),
        "scipy": ("https://docs.scipy.org/doc/scipy/", None),
        "numpy": ("https://numpy.org/doc/stable/", None),
    }
    if packages is None:
        return mapping
    return {k: v for k, v in mapping.items() if k in packages}

