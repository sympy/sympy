from __future__ import annotations
from sympy import Symbol
from sympy.logic.modal import KripkeModel, SemanticEvaluator, Box, Diamond, AgentBox, CommonKnowledge

print("Example 7: Semantic Evaluation and Expressive Operators")

# Kripke Model setup
W = {'w1', 'w2', 'w3'}
R = {'w1': {'w2', 'w3'}}
p = Symbol('p')
V = {
    ('w2', p): True,
    ('w3', p): False
}
model = KripkeModel(W, R, V)
evaluator = SemanticEvaluator(model)

# Evaluation
print(f"Is Diamond(p) true at w1? {evaluator.evaluate(Diamond(p), 'w1')}")
print(f"Is Box(p) true at w1? {evaluator.evaluate(Box(p), 'w1')}")

# Expressive Multi-Agent Operators
agent_box = AgentBox('Alice', p)
common_knowledge = CommonKnowledge('GroupA', p)
print(f"AgentBox modality: {agent_box.modality}")
print(f"CommonKnowledge modality: {common_knowledge.modality}")
