#!/usr/bin/env python

import os

from Cython.Compiler.Main import compile

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext

source_root = os.path.dirname(__file__)

compiled_modules = [
    "sympy.polys.densearith",
    "sympy.polys.densebasic",
    "sympy.polys.densetools",
    "sympy.polys.euclidtools",
    "sympy.polys.factortools",
    "sympy.polys.galoistools",
    "sympy.polys.monomialtools",
    "sympy.polys.orthopolys",
    "sympy.polys.specialpolys",
    "sympy.polys.sqfreetools",
]

extensions = []

for module in compiled_modules:
    source_file = os.path.join(source_root, *module.split('.')) + ".py"

    print("Compiling module %s ..." % module)
    result = compile(source_file)

    if result.c_file is None:
        raise RuntimeError("failed to compile %s" % module)

    extensions.append(
        Extension(module, sources=[str(result.c_file)],
            extra_compile_args=['-O2', '-Wall'],
        )
    )

setup(
    name        = "SymPy",
    packages    = [
        "sympy",
        "sympy.polys",
    ],
    cmdclass    = {
        "build_ext": build_ext
    },
    ext_modules = extensions
)

