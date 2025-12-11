"""
Benchmark for DomainMatrix.to_DM(extension=True) performance.
This script demonstrates the massive performance gap between pure Python
and python-flint when computing primitive elements for high-degree extensions.
"""
import sys
import time
import cProfile
import pstats
from io import StringIO
from sympy import Matrix, sqrt
from sympy.polys.domains import QQ
# --- CONFIGURATION ---
# N=6 creates a field of degree 2^6 = 64.
# In pure Python, this usually takes 10-30 seconds.
# With FLINT, this takes < 0.1 seconds.
# Warning: Setting N >= 12 without FLINT may freeze your computer for minutes.
NUM_SURDS = 10

def get_matrix(n):
    """Generates a matrix with n independent square roots."""
    # Generators: sqrt(2), sqrt(3), sqrt(5)...
    surds = [sqrt(i) for i in range(2, 2 + n)]
    # Create a simple column vector matrix to force extension computation
    return Matrix([[s] for s in surds])

def run_benchmark():
    print(f"--- BENCHMARK CONFIGURATION ---")
    print(f"Generators (N): {NUM_SURDS}")
    print(f"Field Degree:   2^{NUM_SURDS} = {2**NUM_SURDS}")
    # Check if FLINT is currently active in SymPy
    # SymPy's QQ domain will have a library attribute if FLINT is used
    try:
        from sympy.polys.domains.qq import QQ_flint
        has_flint = isinstance(QQ, type(QQ_flint)) or 'flint' in str(type(QQ)).lower()
    except ImportError:
        has_flint = False
    # Explicit check: libraries like python-flint might be installed but not loaded
    # by SymPy unless specifically triggered or configured.
    if 'flint' in sys.modules:
        status = "INSTALLED & ACTIVE (Fast Path)"
    else:
        status = "NOT INSTALLED (Slow Pure-Python Path)"
    print(f"Backend:        {status}")
    print("-" * 60)
    mat = get_matrix(NUM_SURDS)
    print(f"Running to_DM(extension=True)...")
    start_time = time.time()
    # --- CRITICAL SECTION ---
    dm = mat.to_DM(extension=True)
    # ------------------------
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Done.")
    print(f"Elapsed Time:   {elapsed:.4f} seconds")
    print(f"Result Domain:  {dm.domain}")
    print("-" * 60)

def run_profiler():
    """
    Runs cProfile to determine the exact CAUSE of the slowness.
    """
    print("\n--- PROFILING ANALYSIS (DETERMINING THE CAUSE) ---")
    print("Capturing internal function calls...")
    mat = get_matrix(NUM_SURDS)
    pr = cProfile.Profile()
    pr.enable()
    # Run the operation
    mat.to_DM(extension=True)
    pr.disable()
    s = StringIO()
    # Sort by 'cumulative' time to see which function held the program hostage
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(20) # Print top 20 lines
    print(s.getvalue())

if __name__ == "__main__":
    # 1. Run the simple timing benchmark
    run_benchmark()
    # 2. Run the detailed profiler
    run_profiler()
