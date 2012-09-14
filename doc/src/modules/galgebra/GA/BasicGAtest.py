        a,b,c,d,e = MV.setup('a b c d e')
        MV.set_str_format(1)

        print 'e|(a^b) =',e|(a^b)
        print 'e|(a^b^c) =',e|(a^b^c)
        print 'a*(b^c)-b*(a^c)+c*(a^b) =',a*(b^c)-b*(a^c)+c*(a^b)
        print 'e|(a^b^c^d) =',e|(a^b^c^d)
        print -d*(a^b^c)+c*(a^b^d)-b*(a^c^d)+a*(b^c^d)

        print (a^b)|(c^d)

e|(a^b) = {-(b.e)}a
+{(a.e)}b

e|(a^b^c) = {(c.e)}a^b
+{-(b.e)}a^c
+{(a.e)}b^c

a*(b^c)-b*(a^c)+c*(a^b) = {3}a^b^c

e|(a^b^c^d) = {-(d.e)}a^b^c
+{(c.e)}a^b^d
+{-(b.e)}a^c^d
+{(a.e)}b^c^d

{4}a^b^c^d

{(a.d)*(b.c) - (a.c)*(b.d)}1
