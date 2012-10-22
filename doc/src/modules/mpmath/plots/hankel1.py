# Hankel function H1_n(x) on the real line for n=0,1,2,3
h0 = lambda x: hankel1(0,x)
h1 = lambda x: hankel1(1,x)
h2 = lambda x: hankel1(2,x)
h3 = lambda x: hankel1(3,x)
plot([h0,h1,h2,h3],[0,6],[-2,1])