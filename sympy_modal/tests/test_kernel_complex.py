import pytest
from sympy import Symbol, Implies, And, Or, Not
from sympy_modal.frames import KripkeFrame
from sympy_modal.kernel import TrustedKernel, ProofTerm, ModusPonens
from sympy_modal.operators import Box
from sympy_modal.context import ProofContext, Strategy
from sympy_modal.errors import ProofFailure

def test_prover_modus_ponens_chain():
    ctx = ProofContext(KripkeFrame.K())
    p, q, r = Symbol('p'), Symbol('q'), Symbol('r')
    ctx.assume(Implies(p, q, evaluate=False))
    ctx.assume(Implies(q, r, evaluate=False))
    ctx.assume(p)
    proof = ctx.prove(r, strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == r

def test_prover_backward_chain():
    ctx = ProofContext(KripkeFrame.K())
    p, q, r = Symbol('p'), Symbol('q'), Symbol('r')
    ctx.assume(Implies(q, r, evaluate=False))
    ctx.assume(Implies(p, q, evaluate=False))
    ctx.assume(p)
    proof = ctx.prove(r, strategy=Strategy.Backward)
    assert proof.is_valid and proof.formula == r

def test_prover_and_introduction():
    ctx = ProofContext(KripkeFrame.K())
    p, q = Symbol('p'), Symbol('q')
    ctx.assume(p)
    ctx.assume(q)
    proof = ctx.prove(And(p, q, evaluate=False), strategy=Strategy.Backward)
    assert proof.is_valid and proof.formula == And(p, q, evaluate=False)

def test_prover_and_elimination():
    ctx = ProofContext(KripkeFrame.K())
    p, q = Symbol('p'), Symbol('q')
    ctx.assume(And(p, q, evaluate=False))
    proof = ctx.prove(p, strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == p

def test_prover_or_introduction():
    ctx = ProofContext(KripkeFrame.K())
    p, q = Symbol('p'), Symbol('q')
    ctx.assume(p)
    proof = ctx.prove(Or(p, q, evaluate=False), strategy=Strategy.Backward)
    assert proof.is_valid and proof.formula == Or(p, q, evaluate=False)

def test_prover_necessitation_rule():
    ctx = ProofContext(KripkeFrame.K())
    p = Symbol('p')
    # Tautology p -> p
    taut = Implies(p, p, evaluate=False)
    proof = ctx.prove(Box(taut), strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == Box(taut)

def test_prover_k_axiom_distribution():
    ctx = ProofContext(KripkeFrame.K())
    p, q = Symbol('p'), Symbol('q')
    # Prove Box(p -> q) -> (Box(p) -> Box(q))
    # By assuming Box(p -> q) and Box(p)
    ctx.assume(Box(Implies(p, q, evaluate=False)))
    ctx.assume(Box(p))
    proof = ctx.prove(Box(q), strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == Box(q)

def test_prover_t_axiom():
    ctx = ProofContext(KripkeFrame.T())
    p = Symbol('p')
    ctx.assume(Box(p))
    proof = ctx.prove(p, strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == p

def test_prover_four_axiom():
    ctx = ProofContext(KripkeFrame.S4())
    p = Symbol('p')
    ctx.assume(Box(p))
    proof = ctx.prove(Box(Box(p)), strategy=Strategy.ForwardChain)
    assert proof.is_valid and proof.formula == Box(Box(p))

def test_prover_modal_induction_gl():
    ctx = ProofContext(KripkeFrame.GL())
    p = Symbol('p')
    ctx.assume(Box(Implies(Box(p), p, evaluate=False)))
    proof = ctx.prove(Box(p), strategy=Strategy.ModalInduction)
    assert proof.is_valid and proof.formula == Box(p)
