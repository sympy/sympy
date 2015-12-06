""" Implements the function mml2sympy, a parser of Content MathML into a 
	string with the internal representation of the input.
	A sympy object is obtained using the sympify() function on the resulting string.

	mml2sympy depends on sympy and the lxml library.

	Author: Ezequiel Arceo May
	December 2015
"""

from sympy import *
from sympy.printing import mathml
from lxml import etree
x,y,z = symbols('x y z')

# A dictionary to translate tags into sympy objects's names
mml2sympy_dict = {
'root':r"Pow",
'plus':r"Add",
'minus':r'Add',
'times':r'Mul',
'diff':r'Derivative',
'cn':r'Integer',
'power':r'Pow',
'divide':r"Pow",
'ci':r'Symbol',
'int':r'Integral',
'sum':r'Sum',
'sin':r'sin',
'cos':r'cos',
'tan':r'tan',
'cot':r'cot',
'exp':r'exp',
'arcsin':r'asin',
'arcsinh':r'asinh',
'arccos':r'acos',
'arccosh':r'acosh',
'arctan':r'atan',
'arctanh':r'atanh',
'arccot':r'acot',
'arctan':r'atan2',
'ln':r'log',
'eq':r'Equality',
'neq':r'Unequality',
'geq':r'GreaterThan',
'leq':r'LessThan',
'gt':r'StrictGreaterThan',
'lt':r'StrictLessThan',
        }



def mml2sympy(mmltree):
    sympyres = r""
    # #######################
    # ## Non-atomic elements
    if mmltree.tag == "apply":
        # Handle 'minus' tag
        if mmltree[0].tag == 'minus':
            sympyres += r"Add("
            sympyres += mml2sympy(mmltree[1])
            sympyres += r",Mul(Integer(-1),"
            sympyres += mml2sympy(mmltree[2])
            sympyres += r")"
        # Handle 'divide' tag
        elif mmltree[0].tag == 'divide':
            sympyres += r"Mul("
            sympyres += mml2sympy(mmltree[1])
            sympyres += r",Pow("
            sympyres += mml2sympy(mmltree[2])
            sympyres += r",Integer(-1)"
            sympyres += r")"
        # Handle 'root' tag
        elif mmltree[0].tag == 'root':
            sympyres += r"Pow("
            sympyres += mml2sympy(mmltree[1])
            sympyres += r",Rational(1,2)"
            #sympyres += r")"
        # End: handle 'root' tag
        else: 
            sympyres += mml2sympy_dict[mmltree[0].tag]
            sympyres += r"("    
            for branch in list(mmltree)[1:]:
                sympyres += mml2sympy(branch)
                if branch != list(mmltree)[-1]:
                    sympyres += r","
        sympyres += r")"
	# ###################
	# ## Atomic elements
    else:
        content = mmltree.text
        content_type = type(sympify(content))
        # Handle integer content (.text method) in 'cn' tags
        if (content_type == Integer) or (str(content_type)=="<class 'sympy.core.numbers.One'>"):
            sympyres += r"Integer(" + content + r")"
        # Handle float content (.text method) in 'cn' tags
        elif content_type == Float:
            sympyres += r"Float('" + content + r"', prec = 15)"
        # Handle symbol
        else:
            sympyres += r"Symbol('" + content + r"')"
    return(sympyres)



# EXAMPLE
mye = x-y-z -y*z/(x+z*tan(x**2)) + sin(z)+exp(-3/x)*cos(5*x)-1+sqrt(x/(y**2+z**2))
mml=mathml(mye)
tree = etree.XML(mml)
print(etree.tostring(tree,pretty_print=True))
srepr(mye)
back=mml2sympy(tree)
print(back)
print(sympify(back))
print(mye)
