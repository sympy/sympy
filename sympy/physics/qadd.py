from sympy.physics.qoperations import QAssocOp

class QAdd(QAssocOp):
    name = 'QAdd'

    @property
    def eval_to(self):
        if hasattr(self[0], 'eval_to'):
            return self[0].eval_to
        else:
            return self[0].__class__

    @property
    def identity(self):
        return S.Zero
