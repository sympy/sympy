#!/usr/bin/python

import os,sys

print_cnt = 1

def get_fct_str(s):
    if s == '':
        return('','')
    if s[0] == "'" or s[0] == '"':
        istr = s.find(s[0],1)
        if istr == len(s):
            return('','')
        return('',s[istr+2:])

    if ',' in s:
        icomma = s.find(',')
        ilp = s.find('(')
        irp = s.rfind(')')
        if ilp > -1 and irp > -1:
            if icomma > irp:
                return(s[:icomma],s[icomma+1:])
            nbal = 1
            ibal = ilp+1
            n = len(s)
            while ibal < n:
                if s[ibal] == '(':
                    nbal += 1
                if s[ibal] == ')':
                    nbal -= 1
                if nbal == 0:
                    if ibal+1 == n:
                        return(s[:ibal+1],s[ibal+1:])
                    return(s[:ibal+1],s[ibal+2:])
                ibal += 1
            sys.stderr.write('!!!Unbalanced Parenthesis in '+s+'!!!\n')
            exit(1)
        return(s[:icomma],s[icomma+1:])
    return(s,'')

def print_arg_lst(s):
    var_lst = []
    s = s.replace('print ','')
    while s != '':
        (g,s) = get_fct_str(s)
        if g != '':
            var_lst.append(g)
    return(var_lst)

def get_indent(line):
    strp_line = line.strip()
    i = line.find(strp_line[0])
    indent = line[:i]
    return(indent,strp_line)

def adjust_print(line):
    new_line = ''
    (indent,strp_line) = get_indent(line)
    varlst = print_arg_lst(strp_line)
    for var in varlst:
        new_line = indent+'print '+var+'\n'
    return(new_line)

def adjust_Fmt(line):
    new_line = ''
    (indent,strp_line) = get_indent(line)
    iFmt = strp_line.find('.Fmt')
    line = strp_line[:iFmt]
    new_line = indent+'print '+line+'\n'
    return(new_line)

def assert_print(line):
    global print_cnt
    new_line = ''
    (indent,strp_line) = get_indent(line)
    strp_line = strp_line.replace('print ','')
    varlst = print_arg_lst(strp_line)
    for var in varlst:
        new_line += indent+'assert str('+var+') == @'+str(print_cnt)+'\n'
        print_cnt += 1
    return(new_line)

def assert_Fmt(line):
    global print_cnt
    (indent,strp_line) = get_indent(line)
    iFmt = strp_line.find('.Fmt')
    var = strp_line[:iFmt]
    new_line = indent+'assert str('+var+') == @'+str(print_cnt)+'\n'
    print_cnt += 1
    return(new_line)

example_name = sys.argv[1]
#example_name = 'test.py'
test_name = example_name[:-3]+'_test.py'

pgmfile = open(example_name,'r')
#pgmfile = open('./examples/'+example_name,'r')
pgmlst  = pgmfile.readlines()
pgmfile.close()

pgmstr  = ''
teststr = ''
for line in pgmlst:
    if '#' in line:
        line_strp = line.strip()
        if line_strp[0] == '#':
            pass
        else:
            ipound = line.find('#')
            iquote = line.rfind("'")
            if ipound > iquote:
                line = line[:ipound]+'\n'
            else:
                pass
    if 'Eprint' in line:
        pass
    elif 'Format(' in line:
        pass
    elif 'xpdf(' in line:
        pass
    elif 'example_print' in line:
        pass
    elif 'Print_Function' in line:
        pass
    elif 'Get_Program' in line:
        pass
    elif 'print ' in line and 'metric' in line:
        pass
    elif 'def ' in line:
        if line[-4:-1] == '():':
            teststr += line.replace('def ','def test_')
        else:
            teststr += line
        pgmstr   += line

    elif 'print ' in line:
        if line[-1] == "'" or line[-1] == '"':
            pass
        else:
            new_line = adjust_print(line)
            pgmstr += new_line
            new_line = assert_print(line)
            teststr  += new_line
    elif '.Fmt' in line:
        new_line = adjust_Fmt(line)
        pgmstr += new_line
        new_line = assert_Fmt(line)
        teststr  += new_line
    else:
        pgmstr  += line
        teststr += line

print pgmstr
print teststr


datafile = open('data.py','w')
datafile.write(pgmstr)
datafile.close()

os.system('python data.py > data.txt')

datafile = open('data.txt','r')
datalst = datafile.readlines()
datafile.close()
datalst.reverse()

for data in datalst:
    print_cnt -= 1
    teststr = teststr.replace('@'+str(print_cnt),"'"+data.strip()+"'")

ireturn = teststr.rfind('return')
teststr = teststr[:ireturn+6]

testfile = open(test_name,'w')
testfile.write(teststr)
testfile.close()

os.system('rm data.*')
