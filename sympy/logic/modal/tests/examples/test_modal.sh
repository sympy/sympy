#!/bin/bash
pytest sympy/logic/modal/tests/ --cov=sympy.logic.modal --cov-report=term-missing --cov-fail-under=80

echo "Running Example 1:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example1.py
echo ""

echo "Running Example 2:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example2.py
echo ""

echo "Running Example 3:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example3.py
echo ""

echo "Running Example 4:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example4.py
echo ""

echo "Running Example 5:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example5.py
echo ""

echo "Running Example 6:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example6.py
echo ""

echo "Running Example 7:"
PYTHONPATH=. python sympy/logic/modal/tests/examples/example7.py
echo ""
