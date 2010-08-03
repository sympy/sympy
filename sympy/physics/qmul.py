from sympy.physics.qoperations import QAssocOp, _Qrules_

class QMul(QAssocOp):
    name = 'QMul'
    def _validate_QMul(self, other):
        from sympy.physics.quantum import KetBase, BraBase, Operator, InnerProduct
        print 'validating'
        if isinstance(other, (KetBase, BraBase, Operator, InnerProduct)):
            _Qrules_(self.eval_to, other.__class__)
        elif hasattr(other, 'eval_to'):
            _Qrules_(self.eval_to, other.eval_to)
        return

    def _validate_QAdd(self, other):
        from sympy.physics.quantum import KetBase, BraBase, Operator, InnerProduct
        if isinstance(other, (KetBase, BraBase, Operator, InnerProduct)):
            assert other.__class__ == self.eval_to

    @property
    def eval_to(self):
        previous = 0
        for item in self.args:
            if hasattr(item, 'eval_to'):
                item = item.eval_to
            else:
                item = item.__class__
            if previous:
                previous = _Qrules_(previous, item)
            else:
                previous = item
        return previous
                
