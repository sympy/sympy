import abc

class TheorySolver(metaclass=abc.ABCMeta):
    @classmethod
    @abc.abstractmethod
    def from_encoded_cnf(cls, encoded_cnf, *args, **kwargs):
        """
        Construct a theory solver instance from an EncodedCNF,
        returning a tuple (instance, conflict_clauses).
        """
        pass

    @abc.abstractmethod
    def assert_lit(self, literal):
        """
        Assert a literal/constraint into the theory solver.
        Returns None or a conflict explanation.
        """
        pass

    @abc.abstractmethod
    def check(self):
        """
        Check consistency of the theory solver's current assignments.
        Returns (True, model) if consistent or (False, conflict_clause) if inconsistent.
        """
        pass
