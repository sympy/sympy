#! /usr/bin/env bash

if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then

        cd doc
        make clean
        make html

        cd ../../

        git config --global user.email "sympy@googlegroups.com"
        git config --global user.name "SymPy (Travis CI)"

        git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/sympy/sympy_doc.git  gh-pages > /dev/null

        cd gh-pages
        git remote rm origin
        git remote add origin https://${GH_TOKEN}@github.com/sympy/sympy_doc.git
        rm -rf dev/
        cp -R ../sympy/doc/_build/html dev/
        git add -A dev
        ./generate_indexes.py

        git commit -am "Update dev doc after building $TRAVIS_BUILD_NUMBER"
        git push -fq origin gh-pages > /dev/null
fi
