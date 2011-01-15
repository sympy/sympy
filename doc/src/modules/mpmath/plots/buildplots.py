import os.path
import glob

try:
    import psyco
    psyco.full()
except ImportError:
    pass

for f in glob.glob("*.py"):
    if "buildplots" in f or os.path.exists(f[:-3]+".png"):
        continue
    print "Processing", f
    code = open(f).readlines()
    code = ["from mpmath import *; mp.dps=5"] + code
    for i in range(len(code)):
        l = code[i].rstrip()
        if "cplot(" in l:
            l = l[:-1] + (", dpi=45, file='%s.png', verbose=True)" % f[:-3])
            code[i] = l
        elif "splot(" in l:
            l = l[:-1] + (", dpi=45, file='%s.png')" % f[:-3])
            code[i] = l
        elif "plot(" in l:
            l = l[:-1] + (", dpi=45, file='%s.png')" % f[:-3])
            code[i] = l
    code = "\n".join(code)
    exec code
