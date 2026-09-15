from __future__ import annotations
from sympy.logic.modal import FormalisationInterface

fi = FormalisationInterface()
expr_str = "AlethicBox(Symbol('p'))"
formula = fi.formalise(expr_str)

print(f"Parsed formula type: {type(formula).__name__}")

# The interface can infer the correct frame automatically
sig = fi.resolve_modality(expr_str)
print(f"Inferred modality signatures: {sig.operators}")
