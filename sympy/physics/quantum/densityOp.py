from sympy.physics.quantum.qexpr import QExpr, QuantumError
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.core.expr import Expr
from sympy.printing.str import sstr
from sympy.physics.quantum.hilbert import HilbertSpace
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator, OuterProduct
from sympy.physics.quantum.represent import represent
from sympy.functions.elementary.exponential import log
from sympy.physics.quantum.qexpr import _qsympify_sequence
from sympy.core import sympify

class Density(QExpr):
    """
       Generic Density operator object
       The calling convention will be densityOp([state1, prob1],... [staten, probn])
    """

    _label_separator = u' + '
    _op_priority = 11.0
    
    @classmethod
    def _eval_args(cls, args):
        """
            This will check to make sure that all of the pure state hilbert
            spaces contained within it are contained in the same hilbert space?
        """
        for i in xrange(2):
            if isinstance(args[-i-1], list):
                args = args + (1,)
        return sympify(args)
    
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
        for state in self.state_part:
            ret_val = ret_val + (self.lhs*state[0]*Dagger(state[0])*state[1]*self.rhs,)
        return ret_val

    #Should I move the operators into the density at Mul, or on 
    def __mul__(self, other):
        arg_list = list(self.args)
        arg_list[-1] = self.rhs*other
        return Density(*arg_list)

    def __rmul__(self, other):
        arg_list = list(self.args)
        arg_list[-2] = other*self.lhs
        return Density(*arg_list)    
    
    def getstate(self, index):
        #returns the state in place index
        state = self.args[index][0]
        return state
    
    def getprob(self, index):
        #returns the prob in place index
        prob = self.args[index][1]
        return prob  
    
    @property
    def state_part(self):
        return self.args[:-2]
    
    @property
    def lhs(self):
        return self.args[-2]
    
    @property
    def rhs(self):
        return self.args[-1]
        
    def operate_on(self, operators):
        args_list = list(self.args)
        args_list[-1] = self.rhs*Dagger(operators)
        args_list[-2] = operators*self.lhs
        return Density(*args_list)
    
    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, **options):
        bits = options['nqubits']
        state_m = represent(self.getstate(0), nqubits=bits)
        m = state_m*state_m.H*self.getprob(0)

        for i in xrange(len(self.args)-3):
            state_m = represent(self.getstate(i+1), nqubits=bits)
            m += state_m*state_m.H*self.getprob(i+1)
        m = represent(self.lhs, nqubits=bits)*m*represent(self.rhs, nqubits=bits)
        return m 

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------
    def _apply_density(self, **options):
        from sympy.physics.quantum.applyops import apply_operators
        result = []
        for state in self.state_part:
            result.append([apply_operators(self.lhs*state[0], **options), state[1]])
        result.append(1)
        result.append(1)
        return Density(*result)
            

    def entropy(self,nqubits):
        """
            Von Neumann entropy is defined to be sum(lambda*ln(lambda)) 
            where lambda is the eigenvalue of a matrix
        """
        rho = represent(self, nqubits=nqubits)
        from scipy.linalg import eigvals
        return -sum([(eigen*log(eigen)).evalf() for eigen in eigvals(rho).tolist() if eigen != 0])        
