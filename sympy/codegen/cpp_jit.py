from __future__ import annotations

import ctypes
import os
import subprocess
import hashlib
from sympy import Expr, Symbol, cse
from sympy.printing import ccode

from typing import Any

class CPPJITFunction:
    """Wrapper around ctypes CDLL to run JIT-compiled C++ functions."""
    def __init__(self, lib_path: str, num_args: int, num_results: int):
        self._lib_path = lib_path
        self._lib = ctypes.CDLL(lib_path)
        self._num_args = num_args
        self._num_results = num_results

        # Set argtypes and restype for eval_single
        self._lib.eval_single.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double)
        ]
        self._lib.eval_single.restype = None

        # Set argtypes and restype for eval_vectorized
        self._lib.eval_vectorized.argtypes = [
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double)
        ]
        self._lib.eval_vectorized.restype = None

    def __call__(self, *args) -> float | tuple[float, ...] | list[float] | list[tuple[float, ...]] | Any:
        if len(args) == 0:
            raise ValueError("No arguments provided")

        first_arg = args[0]
        
        # Check if numpy is available and first_arg is/contains a numpy array
        try:
            import numpy as np
            is_numpy = isinstance(first_arg, np.ndarray) or (len(args) > 1 and isinstance(args[0], np.ndarray))
        except ImportError:
            is_numpy = False

        if is_numpy:
            if len(args) == 1:
                arr = args[0]
                if arr.ndim != 2 or arr.shape[1] != self._num_args:
                    raise ValueError(f"Expected 2D numpy array with shape (N, {self._num_args})")
                n = arr.shape[0]
                arr = np.ascontiguousarray(arr, dtype=np.float64)
                c_inputs = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                out_arr = np.empty((n, self._num_results), dtype=np.float64)
                c_outputs = out_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                self._lib.eval_vectorized(n, c_inputs, c_outputs)
                return out_arr if self._num_results > 1 else out_arr.squeeze(-1)
            else:
                arr = np.column_stack([np.ascontiguousarray(a, dtype=np.float64) for a in args])
                n = arr.shape[0]
                c_inputs = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                out_arr = np.empty((n, self._num_results), dtype=np.float64)
                c_outputs = out_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                self._lib.eval_vectorized(n, c_inputs, c_outputs)
                return out_arr if self._num_results > 1 else out_arr.squeeze(-1)

        # Check if doing vectorized/batched evaluation
        if isinstance(first_arg, (list, tuple)) or (hasattr(first_arg, '__len__') and not isinstance(first_arg, (str, bytes))):
            if len(args) == 1:
                # Format 1: 2D list of shape (N, num_args)
                rows = list(first_arg)
                n = len(rows)
                flat_inputs = []
                for row in rows:
                    if len(row) != self._num_args:
                        raise ValueError(f"Each input row must have {self._num_args} elements, got {len(row)}")
                    flat_inputs.extend(row)
            else:
                # Format 2: num_args lists/arrays of size N
                if len(args) != self._num_args:
                    raise ValueError(f"Expected {self._num_args} argument arrays, got {len(args)}")
                n = len(first_arg)
                flat_inputs = []
                for i in range(n):
                    for j in range(self._num_args):
                        flat_inputs.append(args[j][i])

            # Allocate ctypes double arrays
            c_inputs = (ctypes.c_double * len(flat_inputs))(*flat_inputs)
            c_outputs = (ctypes.c_double * (n * self._num_results))()

            self._lib.eval_vectorized(n, c_inputs, c_outputs)

            results = []
            for i in range(n):
                row = tuple(c_outputs[i * self._num_results + r] for r in range(self._num_results))
                results.append(row if self._num_results > 1 else row[0])
            return results
        else:
            # Single evaluation
            if len(args) != self._num_args:
                raise ValueError(f"Expected {self._num_args} arguments, got {len(args)}")

            c_inputs = (ctypes.c_double * self._num_args)(*args)
            c_outputs = (ctypes.c_double * self._num_results)()

            self._lib.eval_single(c_inputs, c_outputs)

            res = tuple(c_outputs[r] for r in range(self._num_results))
            return res if self._num_results > 1 else res[0]


def generate_cpp_source(expressions: list[Expr] | Expr, symbols: list[Symbol]) -> str:
    """Generates optimized C++ source code with CSE."""
    if not isinstance(expressions, (list, tuple, set)):
        expressions = [expressions]

    num_args = len(symbols)
    num_results = len(expressions)

    # Run CSE on the expressions
    replacements, reduced_exprs = cse(expressions)

    # C++ Header
    code = "#include <cmath>\n\n"
    code += "extern \"C\" {\n"

    # Shared evaluation logic template
    eval_body = ""
    for idx, sym in enumerate(symbols):
        eval_body += f"    const double {ccode(sym)} = args[{idx}];\n"

    for var, expr in replacements:
        eval_body += f"    const double {ccode(var)} = {ccode(expr)};\n"

    for idx, expr in enumerate(reduced_exprs):
        eval_body += f"    results[{idx}] = {ccode(expr)};\n"

    # Single eval function
    code += "void eval_single(const double* args, double* results) {\n"
    code += eval_body
    code += "}\n\n"

    # Vectorized eval function
    code += "void eval_vectorized(int n, const double* args, double* results) {\n"
    code += "    #pragma omp simd\n"
    code += "    for (int i = 0; i < n; ++i) {\n"
    code += f"        const double* loop_args = args + i * {num_args};\n"
    code += f"        double* loop_results = results + i * {num_results};\n"
    for line in eval_body.splitlines():
        line_replaced = line.replace("args[", "loop_args[").replace("results[", "loop_results[")
        code += f"    {line_replaced}\n"
    code += "    }\n"
    code += "}\n"

    code += "}\n"
    return code


def cpp_jit(expressions: list[Expr] | Expr, symbols: list[Symbol]) -> CPPJITFunction:
    """Compiles the given SymPy expressions to optimized C++ code and JIT-loads it via ctypes."""
    source_code = generate_cpp_source(expressions, symbols)

    # Create the cache directory relative to codegen package
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cache_dir = os.path.join(os.path.dirname(current_dir), "_jit_cache")
    os.makedirs(cache_dir, exist_ok=True)

    h = hashlib.md5(source_code.encode("utf-8")).hexdigest()
    src_path = os.path.join(cache_dir, f"jit_{h}.cpp")
    lib_path = os.path.join(cache_dir, f"jit_{h}.so")

    with open(src_path, "w", encoding="utf-8") as f:
        f.write(source_code)

    # Compile if library does not exist
    if not os.path.exists(lib_path):
        cmd = ["g++", "-O3", "-ffast-math", "-shared", "-fPIC", src_path, "-o", lib_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"C++ JIT Compilation failed:\n{result.stderr}")

    num_args = len(symbols)
    num_results = len(expressions) if isinstance(expressions, (list, tuple, set)) else 1

    return CPPJITFunction(lib_path, num_args, num_results)
