
import sympy

"""
This module assumes you know your way around linear algebra and python.
It will _not: warn you that something you're doing might not be what
you want to do. It'll rather do it and fail gracefull.

"""

def basic_xyz():
    return sympy.symbols("x y z")

def curl(vectorfield,base_symbols=basic_xyz()):
    """
    Takes a vectorfield and the base symbols to be differentiated with.
    If no base is provided, sympy symbols x,y,z are assumed.
    """
    
    x,y,z=base_symbols
    
    val1=sympy.diff(vectorfield[2],y)-sympy.diff(vectorfield[1],z)
    val2=sympy.diff(vectorfield[0],z)-sympy.diff(vectorfield[2],x)
    val3=sympy.diff(vectorfield[1],x)-sympy.diff(vectorfield[0],y)
    
    return Vector(val1,val2,val3)


def divergence(vec,base_symbols=basic_xyz()):
    """returns the sum of the components of the gradient vector"""
    
    gradv=gradient(vec,base_symbols)
    s=0
    for val in gradv:
        s+=val
        
    return s        

def gradient(expr,base_symbols=basic_xyz()):
    """assumes all free symbols in the expression to be relevant variables.
    e.g. 
    (1,1,1) for f(x,y,z)=x+y+z
    (yz,xz,xy) for f(x,y,z)=x*y*z
    """
    x,y,z=base_symbols
    
    n_exprs=[]
    for s in [x,y,z]:
        n_expr=sympy.diff(expr,s)
        n_exprs.append(n_expr)
    return Vector(*n_exprs)

class Matrix:
    """ 
    3x3 matrix, takes vectors as input    
    """
    def __init__(self,v1,v2,v3):
        self.v1=v1
        self.v2=v2
        self.v3=v3
    def __getitem__(self,ind):
        if ind==0:
            return self.v1
        elif ind==1:
            return self.v2
        elif ind==2:
            return self.v3
        else:
            raise ValueError
    
    def __mul__(self,other):
        if isinstance(other,Matrix):
            new_vs=[]
            for v in other:
                rv=self*v
                new_vs.append(rv)
            M=Matrix(*new_vs)
            return M
        if isinstance(other,Vector):
            i=0
            j=0
            new_values=[]
            while i < 3:
                s=0
                while j < 3:
                    
                    val1=self[j][i]
                    val2=other[j]
                    r=val1*val2
                    s+=r
                    j+=1
                
                new_values.append(s)
                
                j=0
                i+=1
                
            V=Vector(*new_values)
            return V

class Vector:
    
    def __init__(self,x,y,z):
        """
        Works basically like a list or a tuple, except it supports
        common vector funtions.
        """
        
        self.x = x 
        self.y = y 
        self.z = z 

    def subs(self,d=None):
        n_exprs=[]
        for expr in self:
            n_expr=expr.subs(d)
            n_exprs.append(n_expr)
            
        return Vector(*n_exprs)
        
    def cross(self,other):
        new_x=self.y*other.z-self.z*other.y
        new_y=self.z*other.x-self.x*other.z
        new_z=self.x*other.y-self.y*self.x
        return Vector(new_x,new_y,new_z)
        
    def __repr__(self):
        s=(self.x,self.y,self.z)
        s=str(s)
        return s
        
    def __getitem__(self,ind):
        if ind==0:
            return self.x
        elif ind==1:
            return self.y
        elif ind==2:
            return self.z
        else:
            raise ValueError
    def __neg__(self):
        return self.__mul__(-1)
    
    def __iter__(self):
        c=0
        while c <3:
            yield self[c]
            c+=1
    
    def __rtruediv__(self,quotient):
        raise TypeError("You can't devide through a vector") 
    
    def __truediv__(self,quotient):
        l=[i/quotient for i in self]        
        return Vector(*l)
    
    def __rmul__(self,other):
        return self*other
    
    def __mul__(self,other):
        """intended for multiplication of the vector with a scalar"""
        if isinstance(other,Vector):
            raise TypeError("use dot product instead")
        
        new_vals=[]
        for val in self:
            new_vals.append(val*other)
            
        return Vector(*new_vals)
        
        
    def __sub__(self,v2):
        return vector3d(self[0]-v2[0],self[1]-v2[1],self[2]-v2[2])
    
    def __add__(self,v2):
        return Vector(self[0]+v2[0],self[1]+v2[1],self[2]+v2[2])
        
    def dot(self,v2):
        v1=self
        return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    
    def dyadic_polynom(self,v2):
        """
        as shown on
        https://en.wikipedia.org/wiki/Dyadic_tensor
        
        returns a polynom
        """
        s=0
        v1=self
        for val1 in v1:
            for val2 in v2:
                s+=val1*val2
        return s
        
    def outer_tensor_product(self,v2):
        """
        tensor product as shown on 
        https://en.wikipedia.org/wiki/Dyadic_tensor
        the function assumes the vectors are already in the appropriate base.
        
        returns a matrix
        """
        
        vectors=[]
        
        for val2 in v2:
            vs=[]
            for val1 in v1:
                vs.append(val1*val2)
                s+=val1*val2
            v=Vector(*vs)
            vectors.append(v)
        m=Matrix(*vectors)
        return m
        
    
    def length(self):
        v=self
        l=(v[0]**2+v[1]**2+v[2]**2)**0.5
        return l

    def normalize(self):
        v=self
        l=self.length()
        return Vector(v[0]/l,v[1]/l,v[2]/l)

def tests():
    v1=Vector(1,0,0)
    v2=Vector(0,1,0)
    v3=Vector(0,0,1)
    v4=Vector(1,2,3)
    
    add_test(v1,v1)
    
    dot_test(v1,v1)
    dot_test(v1,v2)
    
    length_test(v1)
    normalize_test(v1)
    cross_test(v1,v2)
    
    neg_test(v1)
    mul_test(v1)
    
    
    
    quo_test(v1)
    
    #sym_test(v1)
    
    newbasis=Matrix(v1,v3,v2)
    
    matrix_mul_test(newbasis,v4)
    
    x,y,z=sympy.symbols("x y z")
    
    field=1*x+1*y+1*z
    
    gradient_test(field)
    
    divergence_test(field)
    
    v_sym=Vector(x,y,z)
    
    curl_test(v_sym)
    
    subs_test(v_sym,{"x":1})

def subs_test(v1,d):
    v=v1.subs(d)
    print("substest",v)

def divergence_test(vec):
    r=divergence(vec)
    print("divergence",r)
    
def curl_test(vectorfield):
    v=curl(vectorfield)
    print("rotation",v)
    
def gradient_test(field):
    g=gradient(field)
    print("gradient",g)
    
def matrix_mul_test(matrix,v):
    #print(v)
    nv=matrix*v
    print("matrix mul",nv)

def sym_test(v1):
    print("sym",v1[0].free_symbols)

def neg_test(v1):
    v=-v1
    print("neg",v)
    
def quo_test(v1):
    
    r=v1/2
    print("div2",r)
    try:
        r=1/v1
        print("div",r)
    except TypeError:
        pass
    
    

def mul_test(v1):
    v1=5*v1
    print("mul",v1)
    
def add_test(v1,v2):
    v3=v1+v2
    print("add",v3)
    
def dot_test(v1,v2):
    v3=v1.dot(v2)
    print("dot",v3)
    
def cross_test(v1,v2):
    v3=v1.cross(v2)
    print("cross",v3)
    
def length_test(v1):
    l=v1.length()
    print("length",l)
    
def normalize_test(v1):
    nv1=v1.normalize()
    print("norm",nv1)


if __name__=="__main__":
    tests()
    
    

