#!/bin/bash
pytest sympy_modal/tests/ --cov=sympy_modal --cov-report=term-missing --cov-fail-under=80
