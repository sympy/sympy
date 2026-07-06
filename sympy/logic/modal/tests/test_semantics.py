from __future__ import annotations
from sympy.core.symbol import Symbol
from sympy.logic.modal.operators import Box, Diamond
from sympy.logic.modal.semantics import KripkeModel, SemanticEvaluator

def test_semantic_evaluator():
    p = Symbol('p')

    # Model:
    # w1 -> w2
    # w1 -> w3
    # V(w2, p) = True, V(w3, p) = False

    W = {'w1', 'w2', 'w3'}
    R = {'w1': {'w2', 'w3'}}
    V = {
        ('w2', p): True,
        ('w3', p): False
    }

    model = KripkeModel(W, R, V)
    evaluator = SemanticEvaluator(model)

    # Diamond p should be true at w1 (since true at w2)
    assert evaluator.evaluate(Diamond(p), 'w1') == True

    # Box p should be false at w1 (since false at w3)
    assert evaluator.evaluate(Box(p), 'w1') == False


def test_semantic_next_until():
    from sympy.logic.modal.operators import Next, Until
    p = Symbol('p')

    # Model:
    # w1 -> w2 -> w3
    # w1 accesses w2. w2 accesses w3.
    # V(w1, p)=True, V(w2, p)=True, V(w3, q)=True (others False)

    W = {'w1', 'w2', 'w3'}
    R = {'w1': {'w2'}, 'w2': {'w3'}}
    V = {
        ('w1', p): True,
        ('w2', p): True,
        ('w3', q): True
    }

    model = KripkeModel(W, R, V)
    evaluator = SemanticEvaluator(model)

    # Next p at w1 -> true because p is true at w2
    assert evaluator.evaluate(Next(p), 'w1') == True
    # Next q at w1 -> false because q is false at w2
    assert evaluator.evaluate(Next(q), 'w1') == False

    # p Until q at w1 -> true because p holds at w1, w2, and q holds at w3
    assert evaluator.evaluate(Until(p, q), 'w1') == True

    # q Until p at w1 -> true immediately (because p holds at w1... wait no, target is the SECOND argument)
    # Target is q. p holds at w1, w2. q holds at w3.

    # Test cycle failure
    W2 = {'w1', 'w2'}
    R2 = {'w1': {'w2'}, 'w2': {'w1'}} # infinite cycle
    V2 = {('w1', p): True, ('w2', p): True} # q is never true
    m2 = KripkeModel(W2, R2, V2)
    e2 = SemanticEvaluator(m2)
    # p Until q should fail because q is never reached and there's a cycle
    assert e2.evaluate(Until(p, q), 'w1') == False
