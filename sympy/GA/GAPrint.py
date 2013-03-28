#GA_Print.py

import os,sys,StringIO,re,GAdir
from sympy import C,S,Basic,Expr,Symbol,Matrix
from sympy.printing.printer import Printer
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter,accepted_latex_functions
from sympy.printing.conventions import split_super_sub

if GAdir.GA == 'GA':
    from GAsympy import linear_expand
else:
    from sympy.GA.GAsympy import linear_expand

ZERO      = S(0)

SYS_CMD = {'linux2':{'rm':'rm','evince':'evince','null':' > /dev/null','&':'&'},\
           'win32':{'rm':'del','evince':'','null':' > NUL','&':''},\
           'darwin':{'rm':'rm','evince':'open','null':' > /dev/null','&':'&'}}

ColorCode = {
        'black':        '0;30',         'bright gray':  '0;37',
        'blue':         '0;34',         'white':        '1;37',
        'green':        '0;32',         'bright blue':  '1;34',
        'cyan':         '0;36',         'bright green': '1;32',
        'red':          '0;31',         'bright cyan':  '1;36',
        'purple':       '0;35',         'bright red':   '1;31',
        'yellow':       '0;33',         'bright purple':'1;35',
        'dark gray':    '1;30',         'bright yellow':'1;33',
        'normal':       '0'
}

InvColorCode = dict(zip(ColorCode.values(),ColorCode.keys()))

accepted_latex_functions = ['arcsin','arccos','arctan','sin','cos','tan',
                    'theta','beta','alpha','gamma','sinh','cosh','tanh','sqrt',
                    'ln','log','sec','csc','cot','coth','re','im','frac','root',
                    'arg','zeta']

def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None

class enhance_print:

    normal = ''
    base = ''
    fct = ''
    deriv = ''
    bold = ''

    def __init__(self,base=None,fct=None,deriv=None,on=True):
        if on:
            if 'win' in sys.platform:

                if base == None:
                    enhance_print.base = ColorCode['blue']
                else:
                    enhance_print.base = ColorCode[base]
                if fct == None:
                    enhance_print.fct = ColorCode['red']
                else:
                    enhance_print.fct = ColorCode[fct]
                if deriv == None:
                    enhance_print.deriv = ColorCode['cyan']
                else:
                    enhance_print.deriv = ColorCode[deriv]
                enhance_print.normal = '\033[0m'

            else:

                if base == None:
                    enhance_print.base = ColorCode['dark gray']
                else:
                    enhance_print.base = ColorCode[base]
                if fct == None:
                    enhance_print.fct = ColorCode['red']
                else:
                    enhance_print.fct = ColorCode[fct]
                if deriv == None:
                    enhance_print.deriv = ColorCode['cyan']
                else:
                    enhance_print.deriv = ColorCode[deriv]
                enhance_print.normal = '\033[0m'
            print 'Enhanced Printing is on:'
            print 'Base/Blade color is '+InvColorCode[enhance_print.base]
            print 'Function color is '+InvColorCode[enhance_print.fct]
            print 'Derivative color is '+InvColorCode[enhance_print.deriv]+'\n'

            enhance_print.base  = '\033['+enhance_print.base+'m'
            enhance_print.fct   = '\033['+enhance_print.fct+'m'
            enhance_print.deriv = '\033['+enhance_print.deriv+'m'

    @staticmethod
    def enhance_base(s):
        return(enhance_print.base+s+enhance_print.normal)

    @staticmethod
    def enhance_fct(s):
        return(enhance_print.fct+s+enhance_print.normal)

    @staticmethod
    def enhance_deriv(s):
        return(enhance_print.deriv+s+enhance_print.normal)

    @staticmethod
    def strip_base(s):
        new_s = s.replace(enhance_print.base,'')
        new_s = new_s.replace(enhance_print.normal,'')
        return(new_s)

class GA_Printer(StrPrinter):

    function_names = ('acos','acosh','acot','acoth','arg','asin','asinh',\
                      'atan','atan2','atanh','ceiling','conjugate','cos',\
                      'cosh','cot','coth','exp','floor','im','log','re',\
                      'root','sin','sinh','sqrt','sign','tan','tanh')


    def _print_Function(self, expr):
        name = expr.func.__name__
        args = ", ".join([ self._print(arg) for arg in expr.args ])

        if expr.func.nargs is not None:
            if name in GA_Printer.function_names:
                return(expr.func.__name__ + "(%s)"%self.stringify(expr.args, ", "))

        return enhance_print.enhance_fct("%s" % (name,))

    def _print_Derivative(self, expr):
        diff_args = map(self._print,expr.args)
        return(enhance_print.enhance_deriv('D{%s}' % (diff_args[1],))+'%s' % (diff_args[0],))

    def _print_MV(self,expr):
        if expr.obj == ZERO:
            return('0')
        else:
            if expr.print_blades:
                expr.base_to_blade()
            ostr = expr.get_normal_order_str()
            return(ostr)

    def _print_Vector(self,expr):
        if expr.obj == ZERO:
            return('0')
        else:
            ostr = GA_Printer().doprint(expr.obj)
            ostr = ostr.replace(' ','')
            return(ostr)

prog_str = ''
off_mode = False

def Get_Program(off=False):
    global prog_str,off_mode
    prog_file = open(sys.argv[0],'r')
    prog_str = prog_file.read()
    prog_file.close()
    off_mode = off
    return

def Print_Function():
    global prog_str,off_mode
    if off_mode:
        return
    fct_name = str(sys._getframe(1).f_code.co_name)
    ifct = prog_str.find('def '+fct_name)
    iend = prog_str.find('def ',ifct+4)
    tmp_str = prog_str[ifct:iend-1]
    #tmp_str = tmp_str.replace('print_function()','')
    #tmp_str = tmp_str.replace('def '+fct_name+'():\n','')
    fct_name = fct_name.replace('_',' ')
    if GA_LatexPrinter.latex_flg:
        #print '#Code for '+fct_name
        print '##\\begin{lstlisting}[language=Python,showspaces=false,'+\
              'showstringspaces=false,backgroundcolor=\color{gray},frame=single]'
        print tmp_str
        print '##\\end{lstlisting}'
        print '#Code Output:'
    else:
        print '\n'+80*'*'
        #print '\nCode for '+fct_name
        print tmp_str
        print 'Code output:\n'
    return
