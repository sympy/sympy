from fabric.api import env, local, run, sudo, cd, hide, prefix
from fabric.context_managers import shell_env, prefix
from fabric.operations import put, get
from fabric.contrib.files import append, exists
env.use_ssh_config = True

def prepare():
    prepare_apt()
    prepare_userspace()

def prepare_userspace():
    """
    This can be reverted by executing 'remove_userspace'.
    """
    gitrepos()

def prepare_apt():
    sudo("apt-get -qq update")
    sudo("apt-get -y remove libreadline-dev libreadline6-dev libssl-dev libtinfo-dev manpages-dev python-dbus-dev zlib1g-dev")
    sudo("apt-get -y install git python3 make python-virtualenv zip")

def remove_userspace():
    """
    Deletes (!) the SymPy changes. Use with great care.
    """
    run("rm -rf repos")

def gitrepos():
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone https://github.com/sympy/sympy")

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
        # This currently fails with:
        # NameError: global name 'DistutilsFileError' is not defined
        #run("./setup.py bdist_wininst")

def python3_tarball():
    with cd("repos/sympy"):
        run("bin/use2to3")
        with cd("py3ksympy"):
            run("./setup.py clean")
            run("./setup.py sdist")
            # Currently fails:
            #run("./setup.py bdist_wininst")

def build_docs():
    with cd("repos/sympy"):
        run("virtualenv xx")
        run("source xx/bin/activate; pip install sphinx")
        with cd("doc"):
            run("make clean")
            run("source ../xx/bin/activate; make html")
            with cd("_build"):
                run("mv html sympy-docs-html-0.7.0")
                run("zip -9lr sympy-docs-html-0.7.0.zip sympy-docs-html-0.7.0")


def sympy_copy_release_files():
    with cd("repos/sympy"):
        run("mkdir -p /vagrant/release")
        run("cp dist/* /vagrant/release/")
        run("cp doc/_build/sympy-docs-html-0.7.0.zip /vagrant/release")


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
