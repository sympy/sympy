#!/bin/sh

SITE_PACKAGES=`python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"`
SYMPY_PATH=${SITE_PACKAGES}/sympy

if [ -d "$SYMPY_PATH" ]; then
	rm -r $SYMPY_PATH
fi

mkdir -p $SYMPY_PATH

mv /tmp_sympy/* $SYMPY_PATH
rm -r /tmp_sympy
