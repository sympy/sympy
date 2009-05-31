from sympy.logic.boolalg import conjuncts, to_cnf
from sympy.logic.algorithms.dpll import dpll
from sympy.core import Symbol

class KB(object):
    pass

class PropKB(KB):
    "A KB for Propositional Logic.  Inefficient, with no indexing."

    def __init__(self, sentence=None):
        self.clauses = []
        if sentence:
            self.tell(sentence)

    def tell(self, sentence):
        "Add the sentence's clauses to the KB"
        for c in conjuncts(to_cnf(sentence)):
            if not c in self.clauses: self.clauses.append(c)

    def ask(self, query):
        """TODO: examples"""
        if len(self.clauses) == 0: return False
        query_conjuncts = self.clauses[:]
        query_conjuncts.extend(conjuncts(to_cnf(query)))
        s = set()
        for q in query_conjuncts:
            s = s.union(q.atoms(Symbol))
        return bool(dpll(query_conjuncts, list(s), {}))

    def retract(self, sentence):
        "Remove the sentence's clauses from the KB"
        for c in conjuncts(to_cnf(sentence)):
            if c in self.clauses:
                self.clauses.remove(c)