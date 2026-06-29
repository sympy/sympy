from __future__ import annotations

import time
import numpy as np
import sympy as sp
import os
from sympy.codegen.cpp_jit import cpp_jit

def test_n_pendulum_simulation():
    is_ci = "CI" in os.environ
    N = 3 if is_ci else 5
    
    print("\n" + "="*80)
    print(f"CHAOTIC {N}-SEGMENT PENDULUM DYNAMICS ENGINE & SIMULATION")
    print("="*80)
    
    # 1. Symbolic Derivation of Euler-Lagrange equations of motion
    t_start = time.perf_counter()
    
    theta = [sp.Symbol(f"theta_{i}") for i in range(N)]
    omega = [sp.Symbol(f"omega_{i}") for i in range(N)]
    
    # Coordinates of masses
    x = []
    y = []
    for k in range(N):
        x_k = sum(sp.sin(theta[j]) for j in range(k+1))
        y_k = sum(-sp.cos(theta[j]) for j in range(k+1))
        x.append(x_k)
        y.append(y_k)
        
    # Velocities
    vx = []
    vy = []
    for k in range(N):
        vx_k = sum(sp.cos(theta[j]) * omega[j] for j in range(k+1))
        vy_k = sum(sp.sin(theta[j]) * omega[j] for j in range(k+1))
        vx.append(vx_k)
        vy.append(vy_k)
        
    # Kinetic and Potential Energy
    T = 0.5 * sum(vx[k]**2 + vy[k]**2 for k in range(N))
    g = 9.81
    V = g * sum(y[k] for k in range(N))
    
    L = T - V
    
    # Derive Mass Matrix M and Force Vector C: M * alpha + C = 0
    # d/dt (dL/d(omega_i)) - dL/d(theta_i) = 0
    M = sp.Matrix.zeros(N, N)
    C = sp.Matrix.zeros(N, 1)
    
    for i in range(N):
        p_i = sp.diff(L, omega[i])
        for j in range(N):
            M[i, j] = sp.diff(p_i, omega[j])
        C[i, 0] = sum(sp.diff(p_i, theta[j]) * omega[j] for j in range(N)) - sp.diff(L, theta[i])
        
    M = sp.simplify(M)
    C = sp.simplify(C)
    
    t_derive = time.perf_counter() - t_start
    print(f"[*] Symbolic Euler-Lagrange equations derived in: {t_derive:.4f} seconds")
    
    # 2. C++ JIT compilation
    inputs = theta + omega
    outputs = list(M) + list(C)
    
    t_start = time.perf_counter()
    f_jit = cpp_jit(outputs, inputs, force_compile=(not is_ci))
    t_compile = time.perf_counter() - t_start
    print(f"[*] SymPyJIT C++ compilation time: {t_compile:.4f} seconds")
    
    # 3. Compile lambdify as baseline
    from sympy import lambdify
    t_start = time.perf_counter()
    f_lamb = lambdify(inputs, outputs, 'numpy')
    t_lamb_compile = time.perf_counter() - t_start
    print(f"[*] SymPy lambdify (NumPy) compilation time: {t_lamb_compile:.4f} seconds")
    
    # Initial state: small angles, zero velocities
    theta_init = np.array([0.1 * (i + 1) for i in range(N)])
    omega_init = np.zeros(N)
    
    dt = 0.001
    steps = 1000 if is_ci else 150000 # 150,000 physical steps (1,000 in CI)
    
    # RK4 Integration Helper
    def rk4_step(func, t_val, w_val, dt):
        def get_derivatives(cur_t, cur_w):
            state = list(cur_t) + list(cur_w)
            out = func(*state)
            M_flat = out[:N*N]
            C_flat = out[N*N:]
            M_mat = np.array(M_flat).reshape((N, N))
            C_vec = np.array(C_flat)
            alpha_val = np.linalg.solve(M_mat, -C_vec)
            return np.array(cur_w), alpha_val

        k1_w, k1_a = get_derivatives(t_val, w_val)
        k2_w, k2_a = get_derivatives(t_val + 0.5 * dt * k1_w, w_val + 0.5 * dt * k1_a)
        k3_w, k3_a = get_derivatives(t_val + 0.5 * dt * k2_w, w_val + 0.5 * dt * k2_a)
        k4_w, k4_a = get_derivatives(t_val + dt * k3_w, w_val + dt * k3_a)
        
        new_theta = t_val + (dt / 6.0) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)
        new_omega = w_val + (dt / 6.0) * (k1_a + 2 * k2_a + 2 * k3_a + k4_a)
        return new_theta, new_omega

    # --- JIT Simulation ---
    t_start = time.perf_counter()
    theta_jit, omega_jit = theta_init.copy(), omega_init.copy()
    for _ in range(steps):
        theta_jit, omega_jit = rk4_step(f_jit, theta_jit, omega_jit, dt)
    t_sim_jit = time.perf_counter() - t_start
    print(f"[1] SymPyJIT Simulation execution time: {t_sim_jit:.4f} seconds ({steps/t_sim_jit:,.0f} steps/sec)")
    
    # --- Lambdify Simulation ---
    t_start = time.perf_counter()
    theta_lamb, omega_lamb = theta_init.copy(), omega_init.copy()
    for _ in range(steps):
        theta_lamb, omega_lamb = rk4_step(f_lamb, theta_lamb, omega_lamb, dt)
    t_sim_lamb = time.perf_counter() - t_start
    print(f"[2] SymPy lambdify (NumPy) Simulation execution time: {t_sim_lamb:.4f} seconds ({steps/t_sim_lamb:,.0f} steps/sec)")
    
    # Assert correctness
    assert np.allclose(theta_jit, theta_lamb, atol=1e-5)
    
    print("\n" + "-"*85)
    print(f"{'Simulation Engine':<35} | {'Execution Time (s)':<20} | {'Throughput (steps/s)':<22} | {'Speedup':<10}")
    print("-"*85)
    print(f"{'SymPy lambdify (NumPy)':<35} | {t_sim_lamb:<20.4f} | {steps/t_sim_lamb:<22,.0f} | 1.0x")
    print(f"{'SymPyJIT (Native C++ -O3)':<35} | {t_sim_jit:<20.4f} | {steps/t_sim_jit:<22,.0f} | {t_sim_lamb/t_sim_jit:,.1f}x")
    print("="*85 + "\n")
