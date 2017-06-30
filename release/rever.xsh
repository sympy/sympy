# -*-mode: python; flycheck-mode: nil -*-

$XONSH_SHOW_TRACEBACK = True
$RAISE_SUBPROC_ERROR = True

from collections import defaultdict
from collections.abc import Mapping

from rever.activity import activity
from rever.conda import run_in_conda_env

cd ..

$ACTIVITIES = [
    # 'version_bump',
    '_version',
    'mailmap_update',
    'source_tarball',
    'copy_release_files',
    'test_tarball27',
    'test_tarball33',
    'test_tarball34',
    'test_tarball35',
    'test_tarball36',
    'build_docs',

    # 'tag',
]

# Work around https://github.com/ergs/rever/issues/15
@activity
def _version():
    global version
    version = $VERSION

$TAG_PUSH = False

$VERSION_BUMP_PATTERNS = [
    ('sympy/release.py', r'__version__ = ".*"', r'__version__ = "$VERSION"'),
    ]

# ACTIVITIES

@activity
def mailmap_update():
    with run_in_conda_env(['python=3.6', 'mpmath']):
        ./bin/mailmap_update.py

@activity
def source_tarball():
    with run_in_conda_env(['mpmath', 'python=3.6'], 'sympy-release'):
        # Assumes this is run in Docker and git is already clean
        ./setup.py sdist --keep-temp

@activity
def build_docs():
    with run_in_conda_env(['sphinx=1.3.1', 'docutils=0.12', 'numpy', 'mpmath'],
        envname='sympy-release-docs'):

        cd doc
        make clean
        make html-errors
        make man

        cd _build
        mv html @(tarball_format['html-nozip'])
        zip -9lr @(tarball_format['html']) @(tarball_format['html-nozip'])
        cp @(tarball_format['html']) ../../dist/
        cd ..

        make clean
        make latex

        cd _build/latex
        make
        cp @(tarball_format['pdf-orig']) @("../../../dist/{pdf}".format(**tarball_format))
        cd ..


@activity
def copy_release_files():
    ls dist
    cp dist/* /home/release/

@activity
def test_tarball27():
    test_tarball('2.7')

@activity
def test_tarball33():
    test_tarball('3.3')

@activity
def test_tarball34():
    test_tarball('3.4')

@activity
def test_tarball35():
    test_tarball('3.5')

@activity
def test_tarball36():
    test_tarball('3.6')

# HELPER FUNCTIONS

def test_tarball(py_version):
    """
    Test that the tarball can be unpacked and installed, and that sympy
    imports in the install.
    """
    if py_version not in {'2.7', '3.3', '3.4', '3.5', '3.6'}: # TODO: Add win32
        raise ValueError("release must be one of 2.7, 3.3, 3.4, 3.5, or 3.6 not %s" % py_version)


    with run_in_conda_env(['python=%s' % py_version], 'test-install-%s' % py_version):
        cp @('/home/release/{source}'.format(**tarball_format)) @("releasetar.tar".format(**tarball_format))
        tar xvf releasetar.tar

        cd @("/home/{source-orig-notar}".format(**tarball_format))
        python setup.py install
        python -c "import sympy; print(sympy.__version__)"


def get_tarball_name(file):
    """
    Get the name of a tarball

    file should be one of

    source-orig:       The original name of the source tarball
    source-orig-notar: The name of the untarred directory
    source:            The source tarball (after renaming)
    html:              The name of the html zip
    html-nozip:        The name of the html, without ".zip"
    pdf-orig:          The original name of the pdf file
    pdf:               The name of the pdf file (after renaming)
    """
    doctypename = defaultdict(str, {'html': 'zip', 'pdf': 'pdf'})

    if file in {'source-orig', 'source'}:
        name = 'sympy-{version}.tar.gz'
    elif file == 'source-orig-notar':
        name = "sympy-{version}"
    elif file in {'html', 'pdf', 'html-nozip'}:
        name = "sympy-docs-{type}-{version}"
        if file == 'html-nozip':
            # zip files keep the name of the original zipped directory. See
            # https://github.com/sympy/sympy/issues/7087.
            file = 'html'
        else:
            name += ".{extension}"
    elif file == 'pdf-orig':
        name = "sympy-{version}.pdf"
    else:
        raise ValueError(file + " is not a recognized argument")

    ret = name.format(version=version, type=file,
        extension=doctypename[file])
    return ret

tarball_name_types = {
    'source-orig',
    'source-orig-notar',
    'source',
    'html',
    'html-nozip',
    'pdf-orig',
    'pdf',
    }

# Have to make this lazy so that version can be defined.
class _tarball_format(Mapping):
    def __getitem__(self, name):
        return get_tarball_name(name)

    def __iter__(self):
        return iter(tarball_name_types)

    def __len__(self):
        return len(tarball_name_types)

tarball_format = _tarball_format()
