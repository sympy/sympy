from sympy.physics.qoperations import QAssocOp, QRules
from sympy.physics.qoperations import QRules

class QMul(QAssocOp):
    name = 'QMul'    

    @property
    def eval_to(self):
        previous = 0
        for item in self.args:               
            if previous:
                previous = QRules._rules_QMul(previous, item)
            else:
                previous = item
        return previous

    @property
    def identity(self):
        return S.One
                
