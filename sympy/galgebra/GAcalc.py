#!/usr/bin/python
#GAcalc.py

import os.path

#from GAsympy import *

SYSMOD = sys.modules[__name__]

set_main(SYSMOD)

def ExecLine(line):
    global SYSMOD
    if '?' in line:
        return(False)
    if '=' in line:
        """
        Input is a string in form of
            var = expression
        var is extracted from string and input string is
        executed as python expression with exec with result
        placed in globals() dictionary as globals()[var].
        Result is then printed out.
        """
        split_lst = line.split('=')
        lhs = split_lst[0].strip()
        try:
            exec line in globals()
            sys.stdout.write('Out: '+lhs+' = '+str(globals()[lhs])+'\n')
        except:
            print 'Illegal Statement: '+line
            return(False)
    else:
        line = 'dummy = '+line
        line = line.strip()
        if line != 'dummy =':
            try:
                exec line in globals()
                outstr = str(globals()['dummy'])
            except:
                print 'Illegal Statement: '+line
                return(False)
            if outstr != 'None':
                sys.stdout.write('Out: '+outstr+'\n')
    return(True)


if __name__ == '__main__':

    playfile = ''
    savefile = ''
    playfilename = ''

    if len(sys.argv) > 1:
        playfilename = sys.argv[1]
        if os.path.isfile(playfilename):
            playfile = open(playfilename,'r')
        savefile = open('tmptmp','w')

    MV.setup('a0 a1 a2',debug=0)

    line = 'x'

    if playfile != '':
        while line != '':
            line = playfile.readline()
            if line != '':
                ExecLine(line)
                savefile.write(line)

    line = 'x'

    while line != '':
        sys.stdout.write('In: ')
        line = sys.stdin.readline()

        if ExecLine(line):
            if savefile != '':
                savefile.write(line)
        else:
            if savefile != '':
                savefile.close()
                os.system('mv tmptmp '+playfilename)
            break

