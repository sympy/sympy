#Implementation the trigsimp algorithm by fu et al
#Simplifying trigonometric expressions

#The algorithm explanation: http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/DESTIME2006/DES_contribs/Fu/simplification.pdf

#Dimitar Vlahovski @ Technological School "Electronic systems"
#30.11.2011

import re
from sympy import sympify

#-------Transformation rule 0 - simplification of rational polynomials
#Trying to simplify the expression, combine things like 3*x + 2*x, etc.

def TR0(input):
    returnString = sympify(input)
    returnString = simplify(returnString)
    returnString = powsimp(returnString)
    returnString = ratsimp(returnString)
    returnString = combsimp(returnString)

    return str(returnString)

#-------Transformation rule 1 - replace sec, csc with 1/cos, 1/sin 

def TR1(input): 
    returnString = re.sub("sec\((\w*?)\)", "1/cos(\g<1>)", input)
    returnString = re.sub("csc\((\w*?)\)", "1/sin(\g<1>)", returnString)
    
    return returnString
    
#-------Transformation rule 2 - replace tan,cot with sin/cos, cos/sin

def TR2(input):
    returnString = re.sub("tan\((\w*?)\)", "sin(\g<1>)/cos(\g<1>)", input)
    returnString = re.sub("cot\((\w*?)\)", "cos(\g<1>)/sin(\g<1>)", returnString)
 
    return returnString

#-------Transformation rule 3 - induced formula : example sin(-a) = -sin(a)

def TR3(input):
    #Negative argument
    returnString = re.sub("sin\(\-(\w*?)\)", "(-sin(\g<1>))", input)
    returnString = re.sub("cos\(\-(\w*?)\)", "cos(\g<1>)", returnString)
    returnString = re.sub("tan\(\-(\w*?)\)", "(-tan(\g<1>))", returnString)
    returnString = re.sub("cot\(\-(\w*?)\)", "(-cot(\g<1>))", returnString)

    #Argument of type: pi/2 - angle
    returnString = re.sub("sin\(pi/2\-(\w*?)\)", "cos(\g<1>)", returnString) 
    returnString = re.sub("cos\(pi/2\-(\w*?)\)", "sin(\g<1>)", returnString)
    returnString = re.sub("tan\(pi/2\-(\w*?)\)", "cot(\g<1>)", returnString)
    returnString = re.sub("cot\(pi/2\-(\w*?)\)", "tan(\g<1>)", returnString)

    #Argument of type: -angle + pi/2
    returnString = re.sub("sin\(\-(\w*?)\+pi/2\)", "cos(\g<1>)", returnString) 
    returnString = re.sub("cos\(\-(\w*?)\+pi/2\)", "sin(\g<1>)", returnString)
    returnString = re.sub("tan\(\-(\w*?)\+pi/2\)", "cot(\g<1>)", returnString)
    returnString = re.sub("cot\(\-(\w*?)\+pi/2\)", "tan(\g<1>)", returnString)

    #Argument of type: pi/2 + angle
    returnString = re.sub("sin\(pi/2\+(\w*?)\)", "cos(\g<1>)", returnString) 
    returnString = re.sub("cos\(pi/2\+(\w*?)\)", "(-sin(\g<1>))", returnString)
    returnString = re.sub("tan\(pi/2\+(\w*?)\)", "(-cot(\g<1>))", returnString)
    returnString = re.sub("cot\(pi/2\+(\w*?)\)", "(-tan(\g<1>))", returnString)
    
    #Argument of type: angle + pi/2
    returnString = re.sub("sin\((\w*?)\+pi/2\)", "cos(\g<1>)", returnString) 
    returnString = re.sub("cos\((\w*?)\+pi/2\)", "(-sin(\g<1>))", returnString)
    returnString = re.sub("tan\((\w*?)\+pi/2\)", "(-cot(\g<1>))", returnString)
    returnString = re.sub("cot\((\w*?)\+pi/2\)", "(-tan(\g<1>))", returnString)

    #Argument of type: pi - angle
    returnString = re.sub("sin\(pi\-(\w*?)\)", "sin(\g<1>)", returnString)  
    returnString = re.sub("cos\(pi\-(\w*?)\)", "(-cos(\g<1>))", returnString)
    returnString = re.sub("tan\(pi\-(\w*?)\)", "(-tan(\g<1>))", returnString)
    returnString = re.sub("cot\(pi\-(\w*?)\)", "(-cot(\g<1>))", returnString)

    #Argument of type: -angle + pi
    returnString = re.sub("sin\(\-(\w*?)\+pi\)", "sin(\g<1>)", returnString)  
    returnString = re.sub("cos\(\-(\w*?)\+pi\)", "(-cos(\g<1>))", returnString)
    returnString = re.sub("tan\(\-(\w*?)\+pi\)", "(-tan(\g<1>))", returnString)
    returnString = re.sub("cot\(\-(\w*?)\+pi\)", "(-cot(\g<1>))", returnString)

    #Argument of type: pi + angle
    returnString = re.sub("sin\(pi\+(\w*?)\)", "(-sin(\g<1>))", returnString) 
    returnString = re.sub("cos\(pi\+(\w*?)\)", "(-cos(\g<1>))", returnString)
    returnString = re.sub("tan\(pi\+(\w*?)\)", "tan(\g<1>)", returnString)
    returnString = re.sub("cot\(pi\+(\w*?)\)", "cot(\g<1>)", returnString)

    #Argument of type: angle + pi
    returnString = re.sub("sin\((\w*?)\+pi\)", "(-sin(\g<1>))", returnString) 
    returnString = re.sub("cos\((\w*?)\+pi\)", "(-cos(\g<1>))", returnString)
    returnString = re.sub("tan\((\w*?)\+pi\)", "tan(\g<1>)", returnString)
    returnString = re.sub("cot\((\w*?)\+pi\)", "cot(\g<1>)", returnString)

    #Argument of type : 2k*pi - angle
    #2k*pi means that it matches every even number pi e.g. 2*pi, 4*pi, 6*pi ...
    returnString = re.sub("sin\(\d*?(0|2|4|6|8)+?\*pi\s*\-\s*(\w*?)\)", "(-sin(\g<2>))", returnString) 
    returnString = re.sub("cos\(\d*?(0|2|4|6|8)+?\*pi\s*\-\s*(\w*?)\)", "cos(\g<2>)", returnString)
    returnString = re.sub("tan\(\d*?(0|2|4|6|8)+?\*pi\s*\-\s*(\w*?)\)", "(-tan(\g<2>))", returnString)
    returnString = re.sub("cot\(\d*?(0|2|4|6|8)+?\*pi\s*\-\s*(\w*?)\)", "(-cot(\g<2>))", returnString)

    #Argument of type: -angle + 2k*pi
    returnString = re.sub("sin\(\-(\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "(-sin(\g<1>))", returnString) 
    returnString = re.sub("cos\(\-(\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "cos(\g<1>)", returnString)
    returnString = re.sub("tan\(\-(\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "(-tan(\g<1>))", returnString)
    returnString = re.sub("cot\(\-(\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "(-cot(\g<1>))", returnString)


    #Argument of type: 2k*pi + angle
    returnString = re.sub("sin\(\d*?(0|2|4|6|8)+?\*pi\s*\+\s*(\w*?)\)", "sin(\g<2>)", returnString) 
    returnString = re.sub("cos\(\d*?(0|2|4|6|8)+?\*pi\s*\+\s*(\w*?)\)", "cos(\g<2>)", returnString)
    returnString = re.sub("tan\(\d*?(0|2|4|6|8)+?\*pi\s*\+\s*(\w*?)\)", "tan(\g<2>)", returnString)
    returnString = re.sub("cot\(\d*?(0|2|4|6|8)+?\*pi\s*\+\s*(\w*?)\)", "cot(\g<2>)", returnString)

    #Argument of type: angle + 2k*pi
    returnString = re.sub("sin\((\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "sin(\g<1>)", returnString) 
    returnString = re.sub("cos\((\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "cos(\g<1>)", returnString)
    returnString = re.sub("tan\((\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "tan(\g<1>)", returnString)
    returnString = re.sub("cot\((\w*?)\s*\+\s*\d*?(0|2|4|6|8)+?\*pi\)", "cot(\g<1>)", returnString)
   
    return returnString

#-------Transformation rule 4 - values of special angles

def TR4(input):
    returnString = re.sub("sin\(0\)", "0", input)
    returnString = re.sub("cos\(0\)", "1", returnString)
    returnString = re.sub("tan\(0\)", "0", returnString)

    returnString = re.sub("sin\(pi/6\)", "1/2", returnString)
    returnString = re.sub("cos\(pi/6\)", "sqrt(3)/2", returnString)
    returnString = re.sub("tan\(pi/6\)", "sqrt(3)/3", returnString)

    returnString = re.sub("sin\(pi/4\)", "sqrt(2)/2", returnString)
    returnString = re.sub("cos\(pi/4\)", "sqrt(2)/2", returnString)
    returnString = re.sub("tan\(pi/4\)", "1", returnString)

    returnString = re.sub("sin\(pi/3\)", "sqrt(3)/2", returnString)
    returnString = re.sub("cos\(pi/3\)", "1/2", returnString)
    returnString = re.sub("tan\(pi/3\)", "sqrt(3)", returnString)

    returnString = re.sub("sin\(pi/2\)", "1", returnString)
    returnString = re.sub("cos\(pi/2\)", "0", returnString)
   
    return returnString

#-------Transformation rule 5 - substitution of sin square

def TR5(input):
    returnString = re.sub("sin\((\w*?)\)\*\*2", "(1 - cos(\g<1>)**2)", input)
    returnString = re.sub("sin\((\w*?)\)\*\*4", "(1 - cos(\g<1>)**2)**2", returnString)
    
    return returnString

#-------Transformation rule 6 - substitution of cos square

def TR6(input):
    returnString = re.sub("cos\((.*?)\)\*\*2", "(1 - sin(\g<1>)**2)", input)
    returnString = re.sub("cos\((.*?)\)\*\*4", "(1 - sin(\g<1>)**2)**2", returnString)

    return returnString

#-------Transformation rule 7 - Lowering the degree of cos square:

def TR7(input):
    returnString = re.sub("cos\((\w*?)\)\*\*2", "(1+cos(2*\g<1>))/2", input)
    
    return returnString
    
#-------Transformation rule 8 - Converting product to sum or difference

def TR8(input):
    returnString = re.sub("sin\((\w*?)\)\*cos\((\w*?)\)", "(1/2*(sin(\g<1> + \g<2>) + sin(\g<1> - \g<2>)))", input) 
    returnString = re.sub("cos\((\w*?)\)\*sin\((\w*?)\)", "(1/2*(sin(\g<1> + \g<2>) - sin(\g<1> - \g<2>)))", returnString) 
    returnString = re.sub("cos\((\w*?)\)\*cos\((\w*?)\)", "(1/2*(cos(\g<1> + \g<2>) + cos(\g<1> - \g<2>)))", returnString)     
    returnString = re.sub("sin\((\w*?)\)\*sin\((\w*?)\)", "(-1/2*(cos(\g<1> + \g<2>) - cos(\g<1> - \g<2>)))", returnString) 
    
    return returnString

#-------Transformation rule 9 - Converting sum or difference to product

def TR9(input):
    returnString = re.sub("sin\((\w*?)\)\s*\+\s*sin\((\w*?)\)", "2*sin((\g<1> + \g<2>)/2) * cos((\g<1> - \g<2>)/2)", input)  
    returnString = re.sub("sin\((\w*?)\)\s*\-\s*sin\((\w*?)\)", "2*cos((\g<1> + \g<2>)/2) * sin((\g<1> - \g<2>)/2)", returnString)
    returnString = re.sub("cos\((\w*?)\)\s*\+\s*cos\((\w*?)\)", "2*cos((\g<1> + \g<2>)/2) * cos((\g<1> - \g<2>)/2)", returnString)
    returnString = re.sub("cos\((\w*?)\)\s*\-\s*cos\((\w*?)\)", "-2*sin((\g<1> + \g<2>)/2) * sin((\g<1> - \g<2>)/2)", returnString)

    return returnString

#-------Transformation rule 10 - Sum or difference of angles (sin(a+b) = ...)

def TR10(input):
    returnString = re.sub("sin\(\s*(\w*?)\s*\+\s*(\w*?)\s*\)", "(sin(\g<1>)*cos(\g<2>) + cos(\g<1>)*sin(\g<2>))", input)  
    returnString = re.sub("sin\(\s*(\w*?)\s*\-\s*(\w*?)\s*\)", "(sin(\g<1>)*cos(\g<2>) - cos(\g<1>)*sin(\g<2>))", returnString) 
    returnString = re.sub("cos\(\s*(\w*?)\s*\+\s*(\w*?)\s*\)", "(cos(\g<1>)*cos(\g<2>) - sin(\g<1>)*sin(\g<2>))", returnString)
    returnString = re.sub("cos\(\s*(\w*?)\s*\-\s*(\w*?)\s*\)", "(cos(\g<1>)*cos(\g<2>) + sin(\g<1>)*sin(\g<2>))", returnString)
   

    return returnString

#-------Transformation rule 10^-1 - Inverse of Sum or difference of angles:

def TR10_inv(input):
    returnString = re.sub("sin\((?P<arg1>\w*?)\)\s*\*\s*cos\((?P<arg2>\w*?)\)\s*\+\s*cos\((?P=arg1)\)\s*\*\s*sin\((?P=arg2)\)", "sin(\g<1> + \g<2>)", input)  
    returnString = re.sub("sin\((?P<arg1>\w*?)\)\s*\*\s*cos\((?P<arg2>\w*?)\)\s*\-\s*cos\((?P=arg1)\)\s*\*\s*sin\((?P=arg2)\)", "sin(\g<1> - \g<2>)", returnString)  
    returnString = re.sub("cos\((?P<arg1>\w*?)\)\s*\*\s*cos\((?P<arg2>\w*?)\)\s*\-\s*sin\((?P=arg1)\)\s*\*\s*sin\((?P=arg2)\)", "cos(\g<1> + \g<2>)", returnString)  
    returnString = re.sub("cos\((?P<arg1>\w*?)\)\s*\*\s*cos\((?P<arg2>\w*?)\)\s*\+\s*sin\((?P=arg1)\)\s*\*\s*sin\((?P=arg2)\)", "cos(\g<1> - \g<2>)", returnString)  
    return returnString

#-------Transformation rule 11 - Double angle formulas:

def TR11(input):
    returnString = re.sub("sin\((\d*?(0|2|4|6|8)+?)\*(\w*?)\)", "2*sin(\g<1>/2*\g<3>)*cos(\g<1>/2*\g<3>)", input)
    returnString = re.sub("cos\((\d*?(0|2|4|6|8)+?)\*(\w*?)\)", "(cos(\g<1>/2*\g<3>)**2 - sin(\g<1>/2*\g<3>)**2)", returnString)

    return returnString

#-------Transformation rule 12 ---- Sum or difference of tan:

def TR12(input):
    returnString = re.sub("tan\(\s*(\w*?)\s*\+\s*(\w*?)\s*\)", "(tan(\g<1>) + tan(\g<2>)) / (1 - tan(\g<1>)*tan(\g<2>))", input)
    returnString = re.sub("tan\(\s*(\w*?)\s*\-\s*(\w*?)\s*\)", "(tan(\g<1>) - tan(\g<2>)) / (1 + tan(\g<1>)*tan(\g<2>))", returnString)

    return returnString

#-------Transformation rule 13 ---- Product of tan or cot:

def TR13(input):
    returnString = re.sub("tan\((\w*?)\)\*tan\((\w*?)\)", "1 - (tan(\g<1>) + tan(\g<2>))*cot(\g<1> + \g<2>)", input)
    returnString = re.sub("cot\((\w*?)\)\*cot\((\w*?)\)", "1 + (cot(\g<1>) + cot(\g<2>))*cot(\g<1> + \g<2>)", returnString) 

    return returnString

#---------Length of trigonometric expression (count occurances of 'sin' , 'cos', etc.)

def L(input): 
    count = 0
    count += input.count('sin') 
    count += input.count('cos') 
    count += input.count('tan') 
    count += input.count('cot') 
    count += input.count('sec') 
    count += input.count('csc')
    
    return count

#--------Combination transformation rule 1 : ---------------------

def CTR1(input):
    returnString1 = TR5(input)
    returnString1 = TR0(returnString1)

    returnString2 = TR6(input)
    
    returnString2 = TR0(returnString2)
   

    if( (L(returnString1) < L(input)) and (L(returnString1) <= L(returnString2)) ):
       return returnString1
    elif( (L(returnString2) < L(input)) and (L(returnString2) <= L(returnString1)) ):
       return returnString2
    else:
        return input

#--------Combination transformation rule 2: ---------------------

def CTR2(input):
    returnString1 = TR11(input)
    returnString1 = TR5(returnString1)

    returnString2 = TR11(input)
    returnString2 = TR6(returnString2)

    returnString3 = TR11(input)

    returnString1 = TR0(returnString1)
    returnString2 = TR0(returnString2)
    returnString3 = TR0(returnString3)  
   
    if( (L(returnString1) < L(returnString3)) and (L(returnString1) <= L(returnString2)) ):
        return returnString1
    elif( (L(returnString2) < L(returnString3)) and (L(returnString2) <= L(returnString1)) ):
        return returnString2
    else:
        return returnString3

#--------Combination transformation rule 3: ---------------------

def CTR3(input):
    returnString1 = TR8(input)
   

    returnString2 = TR8(input)
    returnString2 = TR10_inv(returnString2)

  

    returnString1 = TR0(returnString1)
    returnString2 = TR0(returnString2)
   
   
    if( (L(returnString2) < L(input)) ):
        return returnString2
    elif( (L(returnString1) < L(input)) ):
        return returnString1
    else:
        return input

#--------Combination transformation rule 4: ---------------------

def CTR4(input):
    returnString1 = TR4(input)
   
    returnString1 = TR0(returnString1)
    
    if( (L(returnString1) < L(input)) ):
        return returnString1
    else:
        return input

#-------------Rule list 1-------------------------

def RL1(input):
    returnString = TR4(input)
    returnString = TR3(returnString)
    returnString = TR4(returnString)
    returnString = TR12(returnString)
    returnString = TR4(returnString)
    returnString = TR13(returnString)
    returnString = TR4(returnString)
    returnString = TR0(returnString)

    return returnString

#-------------Rule list 2-------------------------

def RL2(input):
    returnString = TR4(input)

    returnString = TR3(returnString)
    
    returnString = TR10(returnString)
    
    returnString = TR4(returnString)
    
    returnString = TR3(returnString)
    
    returnString = TR11(returnString)
    
    returnString = TR5(returnString)
    
    returnString = TR7(returnString)
    
    returnString = TR11(returnString)
    
    returnString = TR4(returnString)
    
    returnString = CTR3(returnString)
    
    returnString = TR0(returnString)
    
    returnString = CTR1(returnString)
    
    returnString = TR9(returnString)
    
    returnString = CTR2(returnString)
    
    returnString = TR4(returnString)
    
    returnString = TR9(returnString)
    
    returnString = TR0(returnString)
    
    returnString = TR9(returnString)
    
    returnString = CTR4(returnString)
    
    

    return returnString
   

def my_trigsimp(input):
    """
    ==Usage==
    my_trigsimp(expression) -> reduces expression by using transformation rules given in the fu et al algorithm

    ==Examples===
    No.1
    In [4]: execfile("my_trigsimp.py")

    In [5]: a = sin(50)**2 + cos(50)**2 + sin(pi/6)

    In [6]: my_trigsimp(a)
    Out[6]: 3/2

    No.2
    In [9]: a = sin(100)**4 - cos(50)**2 + sin(50)**2 + 2*cos(100)**2 

    In [10]: my_trigsimp(a)
    Out[10]: 
              2          4     
     2 - 2⋅cos (50) + cos (100)


    """
    input = str(input)
    returnString = TR0(input)
    returnString = TR1(input)
   
    if(returnString.find("tan") or returnString.find("cot")):
        returnString1 = RL1(returnString)
        if(L(returnString1) < L(returnString)):
            returnString = returnString1

    returnString = TR0(returnString)
    returnString = TR2(returnString)
    
    returnString = TR0(returnString)
    

    if(returnString.find("sin") or returnString.find("cos")):
        returnString = RL2(returnString)
        
    
    return sympify(returnString) 
