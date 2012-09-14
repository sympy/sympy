# Airy function Bi(x), Bi'(x) and int_0^x Bi(t) dt on the real line
f = airybi
f_diff = lambda z: airybi(z, derivative=1)
f_int = lambda z: airybi(z, derivative=-1)
plot([f, f_diff, f_int], [-10,2], [-1,2])