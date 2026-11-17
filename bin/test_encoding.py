"""Regression checks for encoding-independent bin scripts."""

import ast
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'bin'


def test_strip_whitespace_utf8(tmp_path):
    source = tmp_path / 'unicode.py'
    source.write_bytes('name = "caf\u00e9"  \r\n'.encode('utf-8'))
    env = {**os.environ, 'PYTHONUTF8': '0', 'PYTHONCOERCECLOCALE': '0',
           'LC_ALL': 'C'}
    subprocess.run([sys.executable, str(BIN / 'strip_whitespace'), str(source)],
                   check=True, env=env, capture_output=True)
    assert source.read_bytes() == 'name = "caf\u00e9"\n'.encode('utf-8')


def test_bin_text_operations_specify_encoding():
    # Keep UTF-8 explicit even on machines where it is the default.
    for filename in ('coverage_doctest.py', 'ask_update.py'):
        tree = ast.parse((BIN / filename).read_text(encoding='utf-8'))
        text_calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                      and isinstance(n.func, ast.Attribute)
                      and n.func.attr in ('read_text', 'write_text')]
        assert text_calls
        assert all(any(k.arg == 'encoding' and isinstance(k.value, ast.Constant)
                       and k.value.value == 'utf-8' for k in n.keywords)
                   for n in text_calls)


def test_executable_check_reads_bytes():
    tree = ast.parse((BIN / 'test_executable.py').read_text(encoding='utf-8'))
    opens = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name) and n.func.id == 'open']
    assert len(opens) == 1
    assert len(opens[0].args) >= 2
    assert isinstance(opens[0].args[1], ast.Constant)
    assert opens[0].args[1].value == 'rb'
