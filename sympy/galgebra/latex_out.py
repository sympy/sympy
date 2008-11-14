#!/usr/bin/python
#latex.py
"""
Convert GAsympy output to LaTeX format
"""

#from string import *
import re

vars = re.compile(r'([A-Za-z]+)([0-9]+)')
subs = re.compile(r'[A-Za-z]+')
xpon = re.compile(r'\*\*([0-9]+)')
greek = re.compile('(alpha|beta|gamma|delta|varepsilon|epsilon|zeta|'+\
         'vartheta|theta|eta|iota|kappa|lambda|mu|nu|xi|varpi|pi|'+\
         'rho|varrho|varsigma|sigma|tau|upsilon|varphi|phi|chi|'+\
         'psi|omega|Gamma|Delta|Theta|Lambda|Xi|Pi|Sigma|Upsilon|'\
         'Phi|Psi|Omega)')

accents = re.compile('([A-Za-z]+)(hat|check|dot|breve|acute|ddot|grave|tilde|mathring|bar|vec)')

GREEK = ['alpha','beta','gamma','delta','varepsilon','epsilon','zeta',\
         'vartheta','theta','eta','iota','kappa','lambda','mu','nu','xi','varpi','pi',\
         'rho','varrho','varsigma','sigma','tau','upsilon','varphi','phi','chi',\
         'psi','omega','Gamma','Delta','Theta','Lambda','Xi','Pi','Sigma','Upsilon',\
         'Phi','Psi','Omega']

ACCENT = ['hat','check','dot','breve','acute','ddot','grave','tilde','mathring',\
          'bar','vec']

def setxpon(x):
    return('^{%s}' % x.group(1))

def setsubscrp(x):
    return('%s_{%s}' % (x.group(1),x.group(2)))

def setaccent(x):
    return('\\%s{%s}' % (x.group(2),x.group(1)))

def setvars(x):
    return('\\%s' % x.group(1))

def LaTeXstr(pstr):
    newpstr = pstr.replace('{','( ')
    newpstr = newpstr.replace('}',' )')
    newpstr = re.sub(vars,setsubscrp,newpstr)
    newpstr = re.sub(accents,setaccent,newpstr)
    newpstr = re.sub(greek,setvars,newpstr)
    newpstr = newpstr.replace('(',r'\lp ')
    newpstr = newpstr.replace(')',r'\rp ')
    newpstr = newpstr.replace('1/2',r'\half ')
    newpstr = newpstr.replace('^',r'\wedge ')
    newpstr = newpstr.replace('|',r'\cdot ')
    newpstr = newpstr.replace('.',r'\cdot ')
    newpstr = re.sub(xpon,setxpon,newpstr)
    newpstr = newpstr.replace('*',' ')
    return(newpstr)

if __name__ == '__main__':

    test = '{(1/2e12^ae45**4+2*psi|c6-1/2ab*alphahat35^c)}'
    print LaTeXstr(test)





