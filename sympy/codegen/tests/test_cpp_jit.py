from __future__ import annotations

import time
import shutil
import pytest

if not shutil.which("g++"):
    pytest.skip("g++ compiler not found, skipping JIT tests.", allow_module_level=True)

np = pytest.importorskip("numpy")

from sympy import symbols, sin, cos, exp, sqrt
from sympy.codegen.cpp_jit import cpp_jit

def test_cpp_jit_basic():
    x, y = symbols('x y')
    # Single expression
    f = cpp_jit(x**2 + y**2, [x, y])
    assert abs(f(3.0, 4.0) - 25.0) < 1e-9

    # Multiple expressions
    f_multi = cpp_jit([x + y, x * y], [x, y])
    assert f_multi(2.0, 3.0) == (5.0, 6.0)


def test_cpp_jit_vectorized():
    x, y = symbols('x y')
    f = cpp_jit(x**2 + y**2, [x, y])

    # Format 1: 2D list of shape (N, num_args)
    res1 = f([[1.0, 2.0], [3.0, 4.0]])
    assert abs(res1[0] - 5.0) < 1e-9
    assert abs(res1[1] - 25.0) < 1e-9

    # Format 2: lists/arrays of size N
    res2 = f([1.0, 3.0], [2.0, 4.0])
    assert abs(res2[0] - 5.0) < 1e-9
    assert abs(res2[1] - 25.0) < 1e-9


def test_cpp_jit_benchmark():
    x, y, z = symbols('x y z')

    # A highly complex multivariate expression with transcendental functions
    expr = sin(x) * cos(y) + exp(z) * sqrt(x**2 + y**2 + z**2 + 1.0) + (x * y) / (z**2 + 1.0)

    # 1,000,000 evaluations
    N = 1000000
    np.random.seed(42)
    x_vals = np.random.rand(N)
    y_vals = np.random.rand(N)
    z_vals = np.random.rand(N)

    # Format inputs for JIT: N x 3 (Python list) and numpy array
    jit_inputs_list = np.column_stack((x_vals, y_vals, z_vals)).tolist()
    jit_inputs_np = np.column_stack((x_vals, y_vals, z_vals))

    print("\n" + "="*80)
    print(f"JIT COMPILER BENCHMARK ENGINE: {N:,} evaluations")
    print("="*80)

    # --- 1. SymPyJIT (C++ Compilation) ---
    t_start = time.perf_counter()
    f_jit = cpp_jit(expr, [x, y, z])
    t_compile = time.perf_counter() - t_start
    print(f"[*] SymPyJIT C++ compilation time: {t_compile:.4f} seconds")

    # --- 2. SymPyJIT with list inputs ---
    t_start = time.perf_counter()
    res_jit_list = f_jit(jit_inputs_list)
    t_jit_list = time.perf_counter() - t_start
    throughput_jit_list = N / t_jit_list
    print(f"[1] SymPyJIT (list) execution time: {t_jit_list:.4f} seconds ({throughput_jit_list:,.0f} evals/sec)")

    # --- 3. SymPyJIT with numpy inputs ---
    t_start = time.perf_counter()
    res_jit_np = f_jit(jit_inputs_np)
    t_jit_np = time.perf_counter() - t_start
    throughput_jit_np = N / t_jit_np
    print(f"[2] SymPyJIT (numpy) execution time: {t_jit_np:.4f} seconds ({throughput_jit_np:,.0f} evals/sec)")

    # --- 4. SymPy lambdify (NumPy) ---
    from sympy import lambdify
    t_start = time.perf_counter()
    f_lamb = lambdify((x, y, z), expr, 'numpy')
    res_lamb = f_lamb(x_vals, y_vals, z_vals)
    t_lamb = time.perf_counter() - t_start
    throughput_lamb = N / t_lamb
    print(f"[3] SymPy lambdify (NumPy) time: {t_lamb:.4f} seconds ({throughput_lamb:,.0f} evals/sec)")

    # --- 5. SymPy evalf (Subset of 100 points, projected to N) ---
    N_evalf = 100
    t_start = time.perf_counter()
    for i in range(N_evalf):
        expr.subs({x: x_vals[i], y: y_vals[i], z: z_vals[i]}).evalf()
    t_evalf_actual = time.perf_counter() - t_start
    t_evalf = t_evalf_actual * (N / N_evalf)
    throughput_evalf = N / t_evalf
    print(f"[4] SymPy evalf (projected) time: {t_evalf:.4f} seconds ({throughput_evalf:,.0f} evals/sec)")

    # Assert correctness
    assert abs(res_jit_list[0] - res_lamb[0]) < 1e-9
    assert abs(res_jit_np[0] - res_lamb[0]) < 1e-9

    # Summary Table
    print("\n" + "-"*85)
    print(f"{'Execution Mode':<35} | {'Time (seconds)':<15} | {'Throughput (evals/s)':<22} | {'Speedup':<10}")
    print("-"*85)
    print(f"{'SymPy evalf (interpreter)':<35} | {t_evalf:<15.4f} | {throughput_evalf:<22,.0f} | 1.0x")
    print(f"{'SymPyJIT (list inputs)':<35} | {t_jit_list:<15.4f} | {throughput_jit_list:<22,.0f} | {t_evalf/t_jit_list:,.1f}x")
    print(f"{'SymPy lambdify (NumPy)':<35} | {t_lamb:<15.4f} | {throughput_lamb:<22,.0f} | {t_evalf/t_lamb:,.1f}x")
    print(f"{'SymPyJIT (Native C++ FFI JIT)':<35} | {t_jit_np:<15.4f} | {throughput_jit_np:<22,.0f} | {t_evalf/t_jit_np:,.1f}x")
    print("="*85 + "\n")
