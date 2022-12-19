from sympy.parsing.mathematica import parse_mathematica
import coverage

cov = coverage.Coverage()
cov.start()

#print(parse_mathematica("-2^-10*x"))
#print(parse_mathematica("-2**x"))
#print(parse_mathematica("Tan[Pi]^Sin[3*Pi/2]"))
#print(parse_mathematica("{[x*(x-2)]-4}+8"))
#print(parse_mathematica("log(x^2)"))
print("parse: d-Log[2*4^-3]")
print("result: ", parse_mathematica("d-log[2*4^-3]"))
#print(parse_mathematica("a^-b*c"))
#print(parse_mathematica("-2^-3"))

cov.stop()
cov.save()
cov.html_report()

#def parse(self, s):
#        s2 = self._from_mathematica_to_tokens(s)
#        s3 = self._from_tokens_to_fullformlist(s2)
#        s4 = self._from_fullformlist_to_sympy(s3)
#        return s4

## zoo == complex infinity
## Testy typu: Stub je objekt, ktorý obsahuje preddefinované údaje a používa ich na odpovedanie na hovory počas testov. Používa sa, keď nemôžeme alebo nechceme zapojiť objekty, ktoré by odpovedali skutočnými údajmi alebo mali nežiaduce vedľajšie účinky.