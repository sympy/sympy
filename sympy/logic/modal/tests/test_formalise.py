from __future__ import annotations
import pytest
from sympy.logic.modal.formalise import FormalisationInterface
from sympy.logic.modal.errors import AmbiguousModalityError, FormalisationError
from sympy.logic.modal.operators import AlethicBox, Box
from sympy.logic.modal.frames import Axiom

def test_formalise_parse():
    fi = FormalisationInterface()
    formula = fi.parse_code("Box(Symbol('p'))")
    assert isinstance(formula, Box)

def test_resolve_modality_explicit():
    fi = FormalisationInterface()
    sig = fi.resolve_modality("AlethicBox(Symbol('p'))")
    assert "alethic" in sig.operators
    assert Axiom.Five in sig.frame.axioms # S5

def test_resolve_modality_ambiguous():
    fi = FormalisationInterface()
    with pytest.raises(AmbiguousModalityError):
         fi.resolve_modality("AlethicBox(EpistemicBox(Symbol('p')))")

def test_resolve_modality_infer_from_axiom():
    fi = FormalisationInterface()
    # T axiom
    sig = fi.resolve_modality("Implies(Box(Symbol('p')), Symbol('p'))")
    # S5 incorporates T
    assert Axiom.T in sig.frame.axioms

def test_resolve_quantifier_scope():
    fi = FormalisationInterface()
    scope = fi.resolve_quantifier_scope("Box(ExistsPredicates(PredicateVariable('P', FunctionType(Universe(0), BoolType())), Symbol('p')))")
    assert not scope.is_ambiguous
    assert len(scope.readings) == 1

def test_resolve_order():
    fi = FormalisationInterface()
    order = fi.resolve_order("ForAllPredicates(PredicateVariable('P', FunctionType(Universe(0), BoolType())), Symbol('p'))")
    assert order.is_second_order
    assert not order.is_first_order

def test_formalise():
    fi = FormalisationInterface()
    res = fi.formalise("AlethicBox(Symbol('p'))")
    assert isinstance(res, AlethicBox)

    err = fi.formalise("NotValidSyntax(*&&^*)")
    assert isinstance(err, FormalisationError)
