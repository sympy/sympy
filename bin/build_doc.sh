#! /usr/bin/env bash

if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then

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
        git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/sympy/sympy_doc.git  gh-pages > /dev/null

        cd gh-pages
        git remote rm origin
        git remote add origin https://${GH_TOKEN}@github.com/sympy/sympy_doc.git > /dev/null
        rm -rf dev/
        cp -R ../sympy/doc/_build/html dev/
        git add -A dev/
        ./generate_indexes.py

        git commit -am "Update dev doc after building $TRAVIS_BUILD_NUMBER"
        echo -e "Pushing commit"
        git push -fq origin gh-pages > /dev/null
fi
