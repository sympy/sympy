from fabric.api import env, local, run, sudo, cd, hide, prefix
from fabric.context_managers import shell_env, prefix
from fabric.operations import put, get
from fabric.contrib.files import append, exists
env.use_ssh_config = True

def prepare():
    prepare_apt()
    checkout_cache()

def prepare_apt():
    sudo("apt-get -qq update")
    sudo("apt-get -y install git python3 make python-virtualenv zip python-dev")
    # Needed to build the docs
    sudo("apt-get -y install graphviz inkscape texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra")

def remove_userspace():
    """
    Deletes (!) the SymPy changes. Use with great care.
    """
    run("rm -rf repos")

def checkout_cache():
    run("rm -rf sympy-cache.git")
    run("git clone --bare https://github.com/sympy/sympy.git sympy-cache.git")

def gitrepos(branch=None):
    if not branch:
        # Use the current branch (of this git repo, not the one in Vagrant)
        branch = local("git rev-parse --abbrev-ref HEAD", capture=True)
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone --reference ../sympy-cache.git https://github.com/sympy/sympy.git")
        if branch != "master":
            with cd("sympy"):
                run("git checkout -t origin/%s" % branch)

def get_sympy_version():
    with cd("repos/sympy"):
        version = run('python -c "import sympy;print sympy.__version__"')
    assert '\n' not in version
    assert ' ' not in version
    assert '\t' not in version
    return version

def test():
    with cd("repos/sympy"):
        run("./setup.py test")

def release(branch=None):
    remove_userspace()
    gitrepos(branch)
    python2_tarball()
    python3_tarball()
    build_docs()
    copy_release_files()

def python2_tarball():
    version = get_sympy_version()
    with cd("repos/sympy"):
        run("git clean -dfx")
        run("./setup.py clean")
        run("./setup.py sdist")
        run("./setup.py bdist_wininst")
        run("mv dist/sympy-{version}.linux-i686.exe dist/sympy-{version}.win32.exe".format(version=version))

def python3_tarball():
    version = get_sympy_version()
    with cd("repos/sympy"):
        run("bin/use2to3")
        with cd("py3k-sympy"):
            run("./setup.py clean")
            run("./setup.py sdist")
            run("mv dist/sympy-{version}.tar.gz dist/sympy-{version}-py3.2.tar.gz".format(version=version))
            run("cp dist/sympy-{version}-py3.2.tar.gz dist/sympy-{version}-py3.3.tar.gz".format(version=version))
            # We didn't test this yet:
            #run("./setup.py bdist_wininst")

def build_docs():
    version = get_sympy_version()
    with cd("repos/sympy"):
        run("mkdir -p dist")
        run("virtualenv xx")
        run("source xx/bin/activate; pip install sphinx==1.1.3 numpy")
        with cd("doc"):
            run("make clean")
            run("source ../xx/bin/activate; make html-errors")
            with cd("_build"):
                run("mv html sympy-docs-html-{version}".format(version=version))
                run("zip -9lr sympy-docs-html-{version}.zip sympy-docs-html-{version}".format(version=version))
                run("cp sympy-docs-html-{version}.zip ../../dist/".format(version=version))
            run("make clean")
            run("source ../xx/bin/activate; make latex")
            with cd("_build"):
                with cd("latex"):
                    run("make")
                    run("cp sympy-{version}.pdf ../../../dist/sympy-docs-pdf-{version}.pdf".format(version=version))

def copy_release_files():
    with cd("repos/sympy"):
        run("mkdir -p /vagrant/release")
        run("cp dist/* /vagrant/release/")
        run("cp py3k-sympy/dist/* /vagrant/release/")

def show_files(file):
    """
    Show the contents of a tarball.

    The current options for file are

    2: The Python 2 tarball
    3: The Python 3 tarball
    2win: The Python 2 Windows installer (Not yet implemented!)
    3win: The Python 3 Windows installer (Not yet implemented!)
    html: The html docs zip
    """
    version = get_sympy_version()
    if file == '2':
        local("tar tf release/sympy-{version}.tar.gz".format(version=version))
    elif file == '3':
        py32 = "sympy-{version}-py3.2.tar.gz".format(version=version)
        py33 = "sympy-{version}-py3.3.tar.gz".format(version=version)
        assert md5(py32) == md5(py33)
        local("tar tf release/" + py32)
    elif file in {'2win', '3win'}:
        raise NotImplementedError("Windows installers")
    elif file == 'html':
        local("unzip -l release/sympy-docs-html-{version}.zip".format(version=version))
    else:
        raise ValueError(file + " is not valid")

def md5(file='*'):
    out = local("md5sum release/" + file, capture=True)
    print out
    return out

# ------------------------------------------------
# Vagrant related configuration

def vagrant():
    vc = _get_vagrant_config()
    # change from the default user to 'vagrant'
    env.user = vc['User']
    # connect to the port-forwarded ssh
    env.hosts = ['%s:%s' % (vc['HostName'], vc['Port'])]
    # use vagrant ssh key
    env.key_filename = vc['IdentityFile'].strip('"')
    # Forward the agent if specified:
    env.forward_agent = vc.get('ForwardAgent', 'no') == 'yes'

def _get_vagrant_config():
    """
    Parses vagrant configuration and returns it as dict of ssh parameters
    and their values
    """
    result = local('vagrant ssh-config', capture=True)
    conf = {}
    for line in iter(result.splitlines()):
        parts = line.split()
        conf[parts[0]] = ' '.join(parts[1:])
    return conf

def restart_network():
    """
    Do this if the VM won't connect to the internet.
    """
    run("sudo /etc/init.d/networking restart")

# ---------------------------------------
# Just a simple testing command:

def uname():
    run('uname -a')
