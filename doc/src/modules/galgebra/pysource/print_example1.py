from printer import Format, xpdf
from ga import Ga
Format()
g3d = Ga('e*x|y|z')
A = g3d.mv('A','mv')
print r'\bm{A} =',A
A.Fmt(2,r'\bm{A}')
A.Fmt(3,r'\bm{A}')
xpdf(paper='letter')
