# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os

def test_CCompilerRunner_environment_variable():
    oldenv = dict(os.environ)
    try:
        os.environ['CC'] = 'ABC'
        from sympy.utilities._compilation.runners import CCompilerRunner
        CCompilerRunner.update_compiler_from_environment_variable()
        assert CCompilerRunner.compiler_dict['default'] == 'ABC'
    finally:
        os.environ.clear()
        os.environ.update(oldenv)

def test_CppCompilerRunner_environment_variable():
    oldenv = dict(os.environ)
    try:
        os.environ['CXX'] = 'ABC'
        from sympy.utilities._compilation.runners import CppCompilerRunner
        CppCompilerRunner.update_compiler_from_environment_variable()
        assert CppCompilerRunner.compiler_dict['default'] == 'ABC'
    finally:
        os.environ.clear()
        os.environ.update(oldenv)

def test_FortranCompilerRunner_environment_variable():
    oldenv = dict(os.environ)
    try:
        os.environ['FC'] = 'ABC'
        from sympy.utilities._compilation.runners import FortranCompilerRunner
        FortranCompilerRunner.update_compiler_from_environment_variable()
        assert FortranCompilerRunner.compiler_dict['default'] == 'ABC'
    finally:
        os.environ.clear()
        os.environ.update(oldenv)
