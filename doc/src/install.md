(installation)=

# Installation

SymPy can be installed on virtually any computer that supports Python.

## From PyPi

The official recommend method of installing Python packages from PyPi is via
pip, with the most basic command being:

```
pip install sympy
```

See the [pip documentation](https://pip.pypa.io/en/stable/index.html) for
variations of this command suitable for your installation needs. Other tools
that pull from PyPi like hatch, poetry, or uv can also be used.

## From anaconda.org

SymPy is packaged for Conda based installers and [available for download on
anaconda.org](https://anaconda.org/search?q=sympy). Install either
[Anaconda](https://www.anaconda.com/products/distribution) or
[Miniconda](https://docs.anaconda.com/miniconda/) and the SymPy distributed
with Anaconda can be installed with:

```
conda install sympy
```

SymPy is also packaged by [Conda Forge](https://conda-forge.org) and if
[Miniforge](https://conda-forge.org/download/) is used, then

```
conda install sympy
```

will install the Conda Forge version of SymPy (which is typically updated
faster than the Anaconda distribution version). You can also install the Conda
Forge version from Anaconda, Miniconda, or Miniforge with:

```
conda install --channel conda-forge sympy
```

Tools such as mamba and pixi can be used to install the SymPy conda package
also.

## From Linux Package Managers

Many Linux distributions package SymPy, for example on Debian based systems
SymPy can be installed with apt:

```
apt install python-sympy
```

or on Fedora based systems, dnf can be used:

```
dnf install sympy
```

## From nightly wheels

We publish a [snapshot of the latest development version of SymPy](
https://anaconda.org/scientific-python-nightly-wheels/sympy) every night as a
pip compatible wheel. You can install the latest version with pip:

```
pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple sympy
```

or with other tools that install wheels.

(installation-git)=
## From Git

If you wish to contribute to SymPy or like to get the latest updates as they
come, install SymPy from git. To download the repository, execute the following
from the command line:

```
git clone https://github.com/sympy/sympy.git
```

To update to the latest version, go into your repository and execute:

```
git pull origin master
```

If you want to install SymPy, but still want to use the git version, you can
run from your repository:

```
python -m pip install --editable .
```

This will cause the installed version to always point to the version in the git
directory.

## Run SymPy

After installation, it is best to verify that your freshly-installed SymPy
works. To do this, start up Python and import the SymPy libraries:

```
$ python
>>> from sympy import *
```

From here, execute some simple SymPy statements like the ones below:

```
>>> x = Symbol('x')
>>> limit(sin(x)/x, x, 0)
1
>>> integrate(1/x, x)
log(x)
```

For a starter guide on using SymPy effectively, refer to the {ref}`intro-tutorial`.

(mpmath-install)=
## mpmath installation

Versions of SymPy prior to 1.0 included [mpmath], but it now depends on it as
an external dependency. If you installed SymPy with pip or conda, it will
already include mpmath. You can manually install mpmath with:

```
pip install mpmath
```

or

```
conda install mpmath
```

to ensure that it is installed.

If you use mpmath via `sympy.mpmath` in your code, you will need to change this
to use just `mpmath`. If you depend on code that does this that you cannot
easily change, you can work around it by doing:

```
import sys
import mpmath
sys.modules['sympy.mpmath'] = mpmath
```

before the code that imports `sympy.mpmath`. It is recommended to change code
that uses `sympy.mpmath` to use `mpmath` directly wherever possible.

## Questions

If you have a question about installation or SymPy in general, feel free to
mail our [mailing list].

If you think there's a bug or you would like to request a feature, please open
an [issue ticket].

[downloads site]: https://github.com/sympy/sympy/releases
[issue ticket]: https://github.com/sympy/sympy/issues
[mailing list]: https://groups.google.com/forum/#!forum/sympy
[mpmath]: https://mpmath.org/
