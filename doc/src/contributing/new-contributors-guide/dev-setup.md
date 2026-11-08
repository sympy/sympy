(devsetup)=

# Setup Development Environment

This guide is intended for people who have never contributed to an open source
project on GitHub before. If you have already completed the steps in this
guide, you do not need to complete them again.

```{note}
This guide is intended for people have never contributed to an open source
project on GitHub before. If you are already familiar with how to contribute
to an open source project on GitHub, go to the [](./workflow-process.md) guide
```

The first step to contributing to the code base is creating your development environment.

```{important}
Each of the steps in this guide only need to be done once. Once you have
completed them, you do not need to repeat them, even if you are making a
second contribution.
```

## Install Git

SymPy is available on [GitHub](https://github.com/sympy/sympy) and uses
[Git](https://git-scm.com/) for source control. The workflow is such that
code is pulled and pushed to and from the main repository. Install the respective version
of Git for your operating system to start development.

**Linux-like systems**:

Install git via your native package management system:

```bash
yum install git
```

or:

```bash
sudo apt-get install git
```

**Windows and macOS**:

The easiest way to get git is to download [GitHub
desktop](https://desktop.github.com/), which will install git, and also
provide a nice GUI (this tutorial will be based on the command line
interface). Note, you may need to go into the GitHub preferences and choose
the "Install Command Line Tools" option to get git installed into the
terminal.

If you do decide to use the GitHub GUI, you should make sure that any "sync
does rebase" option is disabled in the settings.

(configure-git-settings)=
### Configure Your Name and Email in Git

Git tracks who makes each commit by checking the user’s name and email.
In addition, we use this info to associate your commits with your GitHub account.

To set these, enter the code below, replacing the name and email with your own (`--global` is optional).:

```bash
git config --global user.name "Firstname Lastname"
git config --global user.email "your_email@youremail.com"
```

The name should be your actual name, not your GitHub username. Use the email you used for your GitHub account (see [below](dev-setup-create-github-account)).

### (Optional) Configure Git Settings

*This step is not required, but it can make working with git on the command
line easier.*

These global options (i.e. applying to all repositories) are placed in
`~/.gitconfig`. If you want, you can edit this file to enable some handy
shortcuts:

```
[user]
    name = Firstname Lastname
    email = your_email@youremail.com

# Some helpful aliases to save on typing
[alias]
    ci = commit
    di = diff --color-words
    st = status
    co = checkout
    log1 = log --pretty=oneline --abbrev-commit
    logs = log --stat

```

See [Customizing Git - Git Configuration](https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration) for
some more common git configuration options.

## Setup GitHub

Next you will need to setup your GitHub account. Note that all the steps here
only need to be done once. If you already have a GitHub account and have setup
SSH keys, even if it was for a different project than SymPy, you do not need
to do them again.

(dev-setup-create-github-account)=
### Create a GitHub Account

A [GitHub](https://github.com) account is required to contribute to SymPy. If
you have not one yet then [sign up for GitHub](https://github.com/signup). Your
GitHub account is your presence in the open source world, so we recommend
choosing a professional username.

### Setup SSH Keys

To establish a secure connection between your computer and GitHub see detailed
instructions at [Set up Git](https://docs.github.com/get-started/getting-started-with-git/set-up-git) or at [Adding a new SSH key to your GitHub account](https://docs.github.com/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

If you have any problems with SSH access to GitHub, read the troubleshooting
instructions at [Troubleshooting SSH](https://docs.github.com/authentication/troubleshooting-ssh), or
ask us on the [mailing list](https://groups.google.com/g/sympy).

### Fork SymPy

Create your own *fork* of the SymPy project on GitHub. If you have already
done this before, you do not need to do it again.

Go to the [SymPy GitHub repository](https://github.com/sympy/sympy) and click the **Fork** button.

Now you have your own repository for the SymPy project. The address of the
forked project will look something like
`https://github.com/<your-github-username>/sympy`, where
`<your-github-username>` is your GitHub username.

## Get the SymPy Code

It is recommended practice to create a fork of the SymPy project for your development purposes. Create your own fork of the SymPy project (if you have not yet). Go to the SymPy GitHub repository:

```bash
https://github.com/sympy/sympy
```

You will now have a fork at `https://github.com/<your-user-name>/sympy`.

*Note: Replace `<your-github-username>` with your GitHub username.*

Then, on your machine, browse to where you would like to store SymPy, and clone (download) the latest code from SymPy's original repository:

```bash
git clone https://github.com/sympy/sympy
```

For more information about GitHub forking and tuning see: [Fork a Repository](https://docs.github.com/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) and [Pull Requests](https://docs.github.com/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).

### Setup Remote Connection

Remote is a name that git uses locally to store the URL of your forked repository. It makes it easier to collaborate and share your own work without disturbing the original repository. The remote setup can be done in two different ways: via **SSH** or via **HTTPS**.

SSH and HTTPS remote connections are simply two ways of connecting to the forked repository, and are virtually identical to each other. See [Remote Repository](https://docs.github.com/en/get-started/git-basics/about-remote-repositories) for further details.

Move to the cloned sympy repo on your local machine, using the command below:

```
cd sympy
```

### Setup Remote via SSH

After cloning the sympy repository, set up a remote called `github` that points to your forked repository at `git@github.com:<your-github-username>/sympy.git`. To set up the remote using **SSH**, enter the code below:

```bash
git remote add github git@github.com:<your-github-username>/sympy.git
```

*Note: Remote setup via SSH is **only possible** if you have already set up **SSH keys** in your GitHub account. See [Setup SSH Keys](#setup-ssh-keys) to configure them in your GitHub account.*

Then check your remote configuration using the command:

```bash
git remote -v
```

Your setup should be similar to this:

```bash
origin   https://github.com/sympy/sympy (fetch)
origin   https://github.com/sympy/sympy (push)
github git@github.com:<your-github-username>/sympy.git (fetch)
github git@github.com:<your-github-username>/sympy.git (push)
```

### Setup Remote via HTTPS

After cloning the sympy repository, set up a remote called `github` that points to your forked repository at `https://github.com/<your-github-username>/sympy`. To set up the remote using **HTTPS**, enter the code below:

```bash
git remote add github https://github.com/<your-github-username>/sympy.git
```

Then check your remote configuration using the command:

```bash
git remote -v
```

Your setup should be similar to this:

```bash
origin   https://github.com/sympy/sympy (fetch)
origin   https://github.com/sympy/sympy (push)
github https://github.com/<your-github-username>/sympy.git (fetch)
github https://github.com/<your-github-username>/sympy.git (push)
```

## Virtual Environment Setup

You may want to take advantage of using virtual environments to isolate your development version of SymPy from any system wide installed versions, e.g. from `apt-get install python-sympy`.

If you use `conda`, you can use it to create a virtual environment:

```bash
$ conda create -n sympy-dev -c conda-forge --file requirements-dev.txt
```

If you prefer to use `pip` and `venv`, you can use something like

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
```

You can add any other packages to this command that you might find useful for
your contribution, such as the [optional dependencies](../dependencies.md).

You now have a environment that you can use for testing your development copy of SymPy.

Now activate the environment:

```bash
$ conda activate sympy-dev
```
