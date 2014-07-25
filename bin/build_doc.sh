#! /usr/bin/env bash

# This file automatically deploys changes to http://docs.sympy.org/dev/index.html.
# This will happen only when a PR gets merged which is basically when a new commit
# is added to master.
# It requires an access token which should be present in .travis.yml file.
#
# Following is the procedure to get the access token:
#
# $ curl -X POST -u <github_username> -H "Content-Type: application/json" -d\
# "{\"scopes\":[\"public_repo\"],\"note\":\"token for pushing from travis\"}"\
# https://api.github.com/authorizations
#
# It'll give you a JSON response having a key called "token".
#
# $ gem install travis
# $ travis encrypt -r sympy/sympy GH_TOKEN=<token> env.global
#
# This will give you an access token("secure"). This helps in creating an
# environment variable named GH_TOKEN while building.
#
# Add this secure code to .travis.yml as described here http://docs.travis-ci.com/user/encryption-keys/

# Exit on error
set -e

if [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_JOB_NUMBER" == "1" ]; then

        echo "Installing dependencies"
        sudo apt-get install --no-install-recommends graphviz inkscape texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra lmodern librsvg2-bin imagemagick docbook2x
        pip install "sphinx==1.1.3"

        echo -e "Building docs"
        cd doc
        make clean
        make html

        cd ../../
        echo -e "Setting git attributes"
        git config --global user.email "sympy@googlegroups.com"
        git config --global user.name "SymPy (Travis CI)"

        echo -e "Cloning repository"
        git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/sympy/sympy_doc.git  gh-pages > /dev/null 2>&1

        cd gh-pages
        git remote rm origin
        git remote add origin https://${GH_TOKEN}@github.com/sympy/sympy_doc.git > /dev/null 2>&1
        rm -rf dev/
        cp -R ../sympy/doc/_build/html dev/
        git add -A dev/
        ./generate_indexes.py

        git commit -am "Update dev doc after building $TRAVIS_BUILD_NUMBER"
        echo -e "Pushing commit"
        git push -fq origin gh-pages > /dev/null 2>&1
fi
