from sympy.physics.quantum.qexpr import QExpr, QuantumError
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.printing.str import sstr
from sympy.physics.quantum.hilbert import HilbertSpace
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator, OuterProduct
from sympy.physics.quantum.represent import represent

class densityOp(QExpr):
    """
       Generic Density operator object
       The calling convention will be densityOp([state1, prob1],... [staten, probn])
    """

    _label_separator = u' + '
    
    @classmethod
    def eval_args(cls, args):
        """
            This will check to make sure that all of the pure state hilbert
            spaces contained within it are contained in the same hilbert space
            
            generic validation
            
            Possibly mess with the way it is stored
            i.e. |psi><psi| where psi can be anything or |psi1><psi2| where they must be basis states?
        """
        return args
    
    @classmethod
    def _eval_hilbert_space(cls, args):
        """
            The hilbert space is determined by looking at the hilbert space of 
            the internal pure states
            
            for now just return hilbertspace
        """
        return HilbertSpace()
    
    @property    
    def label(self):
        """
            needs to create a tuple of matrix elements in density matrix
            e.g (2*|00><00|, 3*|01><01|)
        
        """
        ret_val = tuple()
        for state in self.args:
            ret_val = ret_val + (state[0][0]*Dagger(state[0][0])*state[1][0],)
        return ret_val

    #Should I move the operators into the density at Mul, or on 
    """def __mul__(self, other):
        if isinstance(other, (Operator, OuterProduct)):
            pass
        else:    
            return Mul(self, other)

    def __rmul__(self, other):
        if isinstance(other, (Operator, OuterProduct)):
            pass
        else:
            return Mul(other, self)
    """
    
    def getstate(self, index):
        #returns the state in place index
        state = self.args[index][0][0]
        return state
    
    def getprob(self, index):
        #returns the prob in place index
        prob = self.args[index][1][0]
        return prob    
    
    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, **options):
        #TODO need to come up with better calling convention
        state_m = represent(self.getstate(0))
        m = state_m*state_m.H*self.getprob(0)

        for i in xrange(len(self.args)-1):
            state_m = represent(self.getstate(i+1))
            m += state_m*state_m.H*self.getprob(i+1)
        return m 

    def entropy(self):
        entropy = 0
        for i in xrange(len(self.args))
            entropy += self.getprob[i]*ln(self.getprob[i])
        return entropy
