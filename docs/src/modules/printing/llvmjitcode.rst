LLVMJITCode Module
==================

This module uses llvmlite to create executable functions from SymPy expressions.
The primary goal is to significantly accelerate the numerical evaluation of
expressions, which is especially useful for functions that need to be called
repeatedly with different arguments (e.g., in numerical integration, optimization, or plotting).

Public API
----------

- `LLVMJitPrinter` : A printer class that converts SymPy expressions into LLVM code.
- `llvm_callable` : Compile a SymPy expression into a callable Python function.

Dependencies
------------

- Requires llvmlite: https://github.com/numba/llvmlite

Example
-------

.. code-block:: python

    from sympy import symbols
    from sympy.printing.llvmjitcode import llvm_callable

    x, y = symbols("x y")
    expr = x**2 + y**2

    f = llvm_callable(expr, [x, y])
    result = f(2, 3)  # returns 13

API Reference
-------------

.. autoclass:: sympy.printing.llvmjitcode.LLVMJitCode
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: sympy.printing.llvmjitcode.LLVMJitPrinter
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: sympy.printing.llvmjitcode.llvm_callable
