**NOTE: The release script is currently in the process of moving from
Vagrant/fabric to Docker/rever. The fabfile.py is left here for reference, but
all release processes should be done with release.sh and rever.xsh.**

# Release

First, make sure that you have done the following things

- Create a release branch. Usually this branch is the same name as the release
  (e.g., "0.7.3"), although no naming convention is enforced on it.

- Change the version in the release branch in sympy/release.py. If you want to
  do a release candidate, change it to a [PEP
  440](https://www.python.org/dev/peps/pep-0440) compliant version like
  0.7.3rc1. Note that setuptools normalizes versions like 0.7.3.rc1 to
  0.7.3rc1, so there will be errors if you do not use the latter form.

- Change the version in master. This way, any additional changes made in master
  will be shown as coming from the right place. The master release should be
  e.g. `0.7.4.dev`, see [PEP 440](https://www.python.org/dev/peps/pep-0440) for
  rules about development version numbers. Note that this version number should
  the next projected version plus the `.dev`.

- Push the release branch up to origin, and make a pull request for it against
  master.

- Create the release notes page for the new release on the wiki. See
  https://github.com/sympy/sympy-bot/issues/26. The easiest way to do this is
  to copy the old release notes to a new page and remove all the changes, and
  update the version number. The formatting on the release notes page is
  important as otherwise the bot will fail, so it is best to do it this way.

It is important to create a new branch because that lets master continue as
normal. The release script will automatically checkout the release branch from
origin, which is why you need to push it (it determines what the release
branch by just looking at what branch you have checked out locally, so make
sure you are on the release branch when you release). It is important to
change the version number because it uses that in naming the tarballs it
creates.

Next, make sure you have Docker installed.

**TODO: Fix the release script to pull sympy/sympy-release from Dockerhub.**

Once you have done these things, execute:

    ./release.sh <BRANCH> <VERSION>

where `<BRANCH>` is the release branch (e.g., `0.7.3`), and `<VERSION>` is the
release version (e.g., `0.7.3rc1`).

On Linux, you may need to use `sudo` to execute this.

This will run all the release scripts. If they are successful, they will
create release tarballs and put them all into a new "release-VERSION"
directory of the current directory. Most likely they will fail the first time,
in which case you will need to investigate why and fix things (e.g., update
authors, run tests, update whitelists in `rever.xsh`, fix setup.py). The whole
script can take about an hour or so to run (depending on how long the tests
take). Every time you re-run the script, it pulls from the branch and runs
everything from scratch.

At the end it will print two things, the list of authors, and the sha256 sums.
Copy the list of authors into the release notes. You should verify that the
sha256 sums of the release files are the same as what are printed.

# Tagging the release

Once you have made the final release files that you plan to upload, be sure
that everything is committed, and that the most recent git HEAD is indeed the
same one that was used to build the files (you can always run the release
script again if you are not sure). Then tag the release with the command

    git tag sympy-VERSION -a

where you should replace `VERSION` with the version (which should be `x.y.z`,
or `x.y.zrcn` for the `n`th release candidate. It is very important to follow
the tag naming conventions.  The `-a` will cause it to prompt for a tag commit
message. Just write something like "SymPy VERSION release".

Then, push up the tag, with

    git push origin sympy-VERSION

Note, once a tag is pushed, that's it. It can't be changed. If you need to
change the tag, you must bump the release number.  So double check that
everything is right before pushing.

# Uploading

**WARNING: This stuff does not fully work yet. Some development on `rever.xsh`
may be required.**

Before you release, you need to push the tag up, as described above.

Release candidates should be uploaded to GitHub only.

    rever VERSION -a GitHub_release

This will create the release on GitHub for the tag, and upload the files to
it.  Do not upload release candidates to PyPI, as `pip` and `easy_install`
will pick them up if you do.

This will prompt you for a username and password the first time you call it.
After that, it will prompt you to generate a token file.  If you don't save
the token to a file, you will need to pass it in as an argument. Releasing is
only supported via OAuth, so using a token is required.

You (obviously) need push access to create a GitHub release.

For final releases, you should upload to both GitHub and PyPI. The command

    rever VERSION -a upload

will do both of these (**TODO: This function has not been translated from the
fabfile yet**).  You will need admin access to the SymPy PyPI project.

Note that if either of these commands fails for some reason, you will very
likely need to go into the web interface and clean some things up before you
can upload again.

# Updating websites

You should now update the websites. Only do this for final releases. The command

    rever VERSION -a update_websites

will update docs.sympy.org and sympy.org (**TODO: This isn't fully translated
from the fabfile yet.**).  You will need to have local clones
of these repos, and push access to them (obviously).  **Note, this command
will commit and push the changes automatically.**

The other website that needs to be updated is SymPy Live. You should make this
as a pull request to the Live repo.

# Updating the Dockerfile

If you change the Dockerfile, you will need to run

    docker build -f Dockerfile . -t sympy/sympy-release

Once you have it working, push the changes up to Dockerhub

    docker push sympy/sympy-release

You'll need access to the sympy org, ask Aaron or Ondřej if you need it.

It is usually not necessary to rebuild the Docker container. The container
first pulls the latest version of the release branch before running rever
(see `pull_and_run_rever.sh`), so unless you modify that script, or change the
packages that are installed in the container, it should not be necessary to
rebuild it.
