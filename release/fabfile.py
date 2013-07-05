from fabric.api import env, local, run, sudo, cd, hide, prefix
from fabric.context_managers import shell_env, prefix
from fabric.operations import put, get
from fabric.contrib.files import append, exists
env.use_ssh_config = True

def prepare(branch="master"):
    prepare_apt()
    prepare_userspace(branch)

def prepare_userspace(branch="master"):
    """
    This can be reverted by executing 'remove_userspace'.
    """
    checkout_cache()
    gitrepos(branch)

def prepare_apt():
    sudo("apt-get -qq update")
    sudo("apt-get -y remove libreadline-dev libreadline6-dev libssl-dev libtinfo-dev manpages-dev python-dbus-dev zlib1g-dev")
    sudo("apt-get -y install git python3 make python-virtualenv zip python-dev")
    # Needed to build the docs
    sudo("apt-get -y install graphviz inkscape texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra")

def remove_userspace():
    """
    Deletes (!) the SymPy changes. Use with great care.
    """
    run("rm -rf repos")

def checkout_cache():
    run("git clone --bare https://github.com/sympy/sympy sympy-cache.git")

def gitrepos(branch="master"):
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone --reference ../sympy-cache.git https://github.com/sympy/sympy")
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

def sympy_test():
    with cd("repos/sympy"):
        run("./setup.py test")

def release():
    python2_tarball()
    python3_tarball()
    build_docs()
    sympy_copy_release_files()

def python2_tarball():
    with cd("repos/sympy"):
        run("git clean -dfx")
        run("./setup.py clean")
        run("./setup.py sdist")
        run("./setup.py bdist_wininst")

def python3_tarball():
    with cd("repos/sympy"):
        run("bin/use2to3")
        with cd("py3k-sympy"):
            run("./setup.py clean")
            run("./setup.py sdist")
            # We didn't test this yet:
            #run("./setup.py bdist_wininst")

def build_docs():
    version = get_sympy_version()
    with cd("repos/sympy"):
        run("virtualenv xx")
        run("source xx/bin/activate; pip install sphinx==1.1.3 numpy")
        with cd("doc"):
            run("make clean")
            run("source ../xx/bin/activate; make html-errors")
            with cd("_build"):
                run("mv html sympy-docs-html-{version}".format(version=version))
                run("zip -9lr sympy-docs-html-{version}.zip sympy-docs-html-{version}".format(version=version))
            run("make clean")
            run("source ../xx/bin/activate; make latex")
            with cd("_build"):
                with cd("latex"):
                    run("make")

def sympy_copy_release_files():
    version=get_sympy_version()
    with cd("repos/sympy"):
        run("mkdir -p /vagrant/release")
        run("cp dist/* /vagrant/release/")
        run("cp doc/_build/sympy-docs-html-{version}.zip /vagrant/release".format(version=version))


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

# ---------------------------------------
# Just a simple testing command:

def uname():
    run('uname -a')
