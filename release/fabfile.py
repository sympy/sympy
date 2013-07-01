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
    hpcmp2_setup()

def prepare_apt():
    sudo("apt-get -qq update")
    sudo("apt-get -y remove libreadline-dev libreadline6-dev libssl-dev libtinfo-dev manpages-dev python-dbus-dev zlib1g-dev")
    sudo("apt-get -y install git make g++ gfortran")

def remove_userspace():
    """
    Deletes (!) the NumPy and Wine changes. Use with great care.
    """
    run("rm -rf repos")

def gitrepos():
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone https://github.com/hashdist/python-hpcmp2")

def hpcmp2_setup():
    with cd("repos/python-hpcmp2"):
        put("config.yml", ".")

def hpcmp2_build():
    with cd("repos/python-hpcmp2"):
        run("./update")

def hpcmp2_check_libs():
    with cd("repos/python-hpcmp2"):
        run("./update --check-libs")

# ------------------------------------------------
# Vagrant related configuration

def vagrant():
    vc = _get_vagrant_config()
    # change from the default user to 'vagrant'
    env.user = vc['User']
    # connect to the port-forwarded ssh
    env.hosts = ['%s:%s' % (vc['HostName'], vc['Port'])]
    # use vagrant ssh key
    env.key_filename = vc['IdentityFile']
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
