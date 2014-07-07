from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga
from sympy.galgebra.mv import Com
Format()

(g4d,a,b,c,d) = Ga.build('a b c d')
print '\\bm{a|(b*c)} =',a|(b*c)
print '\\bm{a|(b^c)} =',a|(b^c)
print '\\bm{a|(b^c^d)} =',a|(b^c^d)
print '\\bm{a|(b^c)+c|(a^b)+b|(c^a)} =',(a|(b^c))+(c|(a^b))+(b|(c^a))
print '\\bm{a*(b^c)-b*(a^c)+c*(a^b)} =',a*(b^c)-b*(a^c)+c*(a^b)
print '\\bm{a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)} =',a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)
print '\\bm{(a^b)|(c^d)} =',(a^b)|(c^d)
print '\\bm{((a^b)|c)|d} =',((a^b)|c)|d
print '\\bm{(a^b)\\times (c^d)} =',Com(a^b,c^d)

xpdf(paper='letter')
