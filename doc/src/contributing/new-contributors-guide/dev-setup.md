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
[Git](https://git-scm.com) for source control. The workflow is such that
code is pulled and pushed to and from the main repository. Install the respective version
of Git for your operating system to start development.

**Linux-like systems**:

Install git via your native package management system:

```
yum install git
```

or:

```
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

### Configure Your Name and Email in Git

Git tracks who makes each commit by checking the user’s name and email.
In addition, we use this info to associate your commits with your GitHub account.

To set these, enter the code below, replacing the name and email with your own (`--global` is optional).:

```
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

See <https://git-scm.com/book/sv/v2/Customizing-Git-Git-Configuration> for
some more common git configuration options.

### (Optional) Tune Your Bash Prompt for Git

*This step is not required, but it can make working with git on the command
line easier, so it is recommended.*

It can be convenient in future to tune the bash prompt to display the current git branch.

The easiest way to do it, is to add the snippet below to your .bashrc or .bash_profile:

```
PS1="[\u@\h \W\$(git branch 2> /dev/null | grep -e '\* ' | sed 's/^..\(.*\)/{\1}/')]\$ "
```

But better is to use `git-completion` from the `git` source. This also has the advantage of adding tab completion to just about every git command. It also includes many other useful features, for example,
promptings. To use `git-completion`, first download the `git` source code (about 27 MiB), then copy
the file to your profile directory:

```
git clone git://git.kernel.org/pub/scm/git/git.git
cp git/contrib/completion/git-completion.bash ~/.git-completion.sh
```

Read instructions in '~/.git-completion.sh'

Note that if you install git from the package manager in many Linux distros, this file is already installed for you.  You can check if it is installed by seeing if tab completion works on git commands (try, e.g., `git commi<TAB>`, or `git log --st<TAB>`). You can also check if the PS1 commands work by doing something like:

```
PS1='\W $(__git_ps1 "%s")\$ '
```

And your command prompt should change to something like:

```
sympy master$
```

Note, it is important to define your PS1 using single quotes ('), not double quotes ("), or else bash will not update the branch name.

## Setup GitHub

Next you will need to setup your GitHub account. Note that all the steps here
only need to be done once. If you already have a GitHub account and have setup
SSH keys, even if it was for a different project than SymPy, you do not need
to do them again.

(dev-setup-create-github-account)=
### Create a GitHub Account

A [GitHub](https://github.com) account is required to contribute to SymPy. If
you have not one yet then sign up at <https://github.com/signup/free>. Your
GitHub account is your presence in the open source world, so we recommend
choosing a professional username.

### Setup SSH Keys

To establish a secure connection between your computer and GitHub see detailed
instructions in <https://help.github.com/en/articles/set-up-git> or at
<https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account>.

If you have any problems with SSH access to GitHub, read the troubleshooting
instructions at <https://help.github.com/en/articles/troubleshooting-ssh>, or
ask us on the [mailing list](https://groups.google.com/group/sympy).

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

Then, on your machine browse to where you would like to store SymPy, and clone (download) the latest code from SymPy's original repository (about 77 MiB):

```bash
$ git clone https://github.com/sympy/sympy
```

Then assign your read-and-write repo to a remote called "github" (replace
`<your-github-username>` with your GitHub username):

```
git remote add github git@github.com:<your-github-username>/sympy.git
```

For more information about GitHub forking and tuning see:
<https://help.github.com/en/articles/about-pull-requests>, <https://help.github.com/en/articles/fork-a-repo>, and <https://help.github.com/en/articles/set-up-git>

After the configuration, your setup should be similar to this:

```bash
$ git remote -v
origin   https://github.com/sympy/sympy (fetch)
origin   https://github.com/sympy/sympy (push)
upstream https://github.com/github/sympy (fetch)
upstream https://github.com/github/sympy (push)
```

## Virtual Environment Setup

You may want to take advantage of using virtual environments to isolate your development version of SymPy from any system wide installed versions, e.g. from `apt-get install python-sympy`.

We recommend using `conda` to create a virtual environment:

```bash
$ conda create -n sympy-dev python=3 mpmath flake8
```

You can add any other packages to this command that you might find useful for
your contribution, such as the [optional dependencies](../dependencies.md).

You now have a environment that you can use for testing your development copy of SymPy.

Now activate the environment:

```bash
$ conda activate sympy-dev
```
