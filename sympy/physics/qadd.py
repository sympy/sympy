from sympy.physics.qoperations import QAssocOp

class QAdd(QAssocOp):
    def _validate_QMul(self, other):
        pass

    def _validate_QAdd(self, other):
        pass

    @property
    def eval_to(self):
        pass
