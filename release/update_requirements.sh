#!/usr/bin/env bash
#
# Run as
#
#   $ ./update_requirements.sh python
#
# where python is the name of the (on PATH) interpreter to use. This will
# create a virtual environment with venv, rerun the pip install command to get
# latest versions meeting the constraints and then use pip list to recreate
# requirements.txt with updated versions.
#
# This will overwrite the requirements.txt file and show the diff

if [ -z $1 ];
then
  python=python3
else
  python=$1
fi

tmp_dir=$(mktemp -d -t venv-XXXX)
venv_dir=$tmp_dir/venv

$python -m venv $venv_dir
. $venv_dir/bin/activate

pip install -U pip wheel

# Installing tensorflow needs a lot of memory. Give at least 2GiB memory for a
# VM or it will fail printing "Collecting tensorflow... Killed"

pip install\
  mpmath\
  'matplotlib>=2.2'\
  'numpy==1.18.5'\
  scipy\
  theano\
  ipython\
  symengine\
  tensorflow\
  cython\
  llvmlite\
  wurlitzer\
  autowrap\
  numexpr\
  'antlr4-python3-runtime==4.7.*'\
  sphinx\
  sphinx-math-dollar\
  #

pip freeze > requirements.txt
git diff requirements.txt
