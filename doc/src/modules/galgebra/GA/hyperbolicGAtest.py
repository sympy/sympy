        print 'Example: non-euclidian distance calculation'

        metric = '0 # #,# 0 #,# # 1'
        X,Y,e = MV.setup('X Y e',metric)
        XdotY = sympy.Symbol('(X.Y)')
        Xdote = sympy.Symbol('(X.e)')
        Ydote = sympy.Symbol('(Y.e)')
        MV.set_str_format(1)
        L = X^Y^e
        B = L*e
        Bsq = (B*B)()
        print 'L = X^Y^e is a non-euclidian line'
        print 'B = L*e =',B
        BeBr =B*e*B.rev()
        print 'B*e*B.rev() =',BeBr
        print 'B^2 =',Bsq
        print 'L^2 =',(L*L)()
        s,c,Binv,M,S,C,alpha = symbols('s c Binv M S C alpha')
        Bhat = Binv*B # Normalize translation generator
        R = c+s*Bhat # Rotor R = exp(alpha*Bhat/2)
        print 's = sinh(alpha/2) and c = cosh(alpha/2)'
        print 'R = exp(alpha*B/(2*|B|)) =',R
        Z = R*X*R.rev()
        Z.expand()
        Z.collect([Binv,s,c,XdotY])
        print 'R*X*R.rev() =',Z
        W = Z|Y
        W.expand()
        W.collect([s*Binv])
        print '(R*X*rev(R)).Y =',W
        M = 1/Bsq
        W.subs(Binv**2,M)
        W.simplify()
        Bmag = sympy.sqrt(XdotY**2-2*XdotY*Xdote*Ydote)
        W.collect([Binv*c*s,XdotY])

        W.subs(2*XdotY**2-4*XdotY*Xdote*Ydote,2/(Binv**2))
        W.subs(2*c*s,S)
        W.subs(c**2,(C+1)/2)
        W.subs(s**2,(C-1)/2)
        W.simplify()
        W.subs(1/Binv,Bmag)
        W = W().expand()
        print '(R*X*R.rev()).Y =',W
        nl = '\n'

        Wd = collect(W,[C,S],exact=True,evaluate=False)
        print 'Wd =',Wd
        Wd_1 = Wd[ONE]
        Wd_C = Wd[C]
        Wd_S = Wd[S]
        print '|B| =',Bmag
        Wd_1 = Wd_1.subs(Bmag,1/Binv)
        Wd_C = Wd_C.subs(Bmag,1/Binv)
        Wd_S = Wd_S.subs(Bmag,1/Binv)
        print 'Wd[ONE] =',Wd_1
        print 'Wd[C] =',Wd_C
        print 'Wd[S] =',Wd_S


        lhs = Wd_1+Wd_C*C
        rhs = -Wd_S*S
        lhs = lhs**2
        rhs = rhs**2
        W = (lhs-rhs).expand()
        W = (W.subs(1/Binv**2,Bmag**2)).expand()
        print 'W =',W
        W = (W.subs(S**2,C**2-1)).expand()
        print 'W =',W
        W = collect(W,[C,C**2],evaluate=False)
        print 'W =',W

        a = W[C**2]
        b = W[C]
        c = W[ONE]

        print 'a =',a
        print 'b =',b
        print 'c =',c

        D = (b**2-4*a*c).expand()
        print 'Setting to 0 and solving for C gives:'
        print 'Descriminant D = b^2-4*a*c =',D
        C = (-b/(2*a)).expand()
        print 'C = cosh(alpha) = -b/(2*a) =',C

Example: non-euclidian distance calculation
L = X^Y^e is a non-euclidian line
B = L*e = X^Y
+{-(Y.e)}X^e
+{(X.e)}Y^e

B*e*B.rev() = {2*(X.Y)*(X.e)*(Y.e) - (X.Y)**2}e

B^2 = -2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2
L^2 = -2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2
s = sinh(alpha/2) and c = cosh(alpha/2)
R = exp(alpha*B/(2*|B|)) = {c}1
+{Binv*s}X^Y
+{-(Y.e)*Binv*s}X^e
+{(X.e)*Binv*s}Y^e

R*X*R.rev() = {Binv*(2*(X.Y)*c*s - 2*(X.e)*(Y.e)*c*s) + Binv**2*((X.Y)**2*s**2
               - 2*(X.Y)*(X.e)*(Y.e)*s**2) + c**2}X
             +{2*Binv*c*s*(X.e)**2}Y
             +{Binv**2*(-2*(X.e)*(X.Y)**2*s**2 + 4*(X.Y)*(Y.e)*(X.e)**2*s**2)
               - 2*(X.Y)*(X.e)*Binv*c*s}e

(R*X*rev(R)).Y = {Binv*s*(-4*(X.Y)*(X.e)*(Y.e)*c + 2*c*(X.Y)**2)
                 + Binv**2*s**2*(-4*(X.e)*(Y.e)*(X.Y)**2 +
                 4*(X.Y)*(X.e)**2*(Y.e)**2 + (X.Y)**3) + (X.Y)*c**2}1

(R*X*R.rev()).Y = S*(-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2)
                  + (X.Y)*Binv*C*(-2*(X.Y)*(X.e)*(Y.e) +
                  (X.Y)**2)**(1/2) + (X.e)*(Y.e)*Binv*(-2*(X.Y)*(X.e)*(Y.e)
                  + (X.Y)**2)**(1/2) -
                  (X.e)*(Y.e)*Binv*C*(-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2)
Wd = {1: (X.e)*(Y.e)*Binv*(-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2),
      S: (-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2),
      C: (X.Y)*Binv*(-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2) -
         (X.e)*(Y.e)*Binv*(-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2)}
|B| = (-2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2)**(1/2)

Wd[ONE] = (X.e)*(Y.e)
Wd[C] = (X.Y) - (X.e)*(Y.e)
Wd[S] = 1/Binv

W = 2*(X.Y)*(X.e)*(Y.e)*C + (X.Y)**2*C**2 + (X.e)**2*(Y.e)**2
    - (X.Y)**2*S**2 + (X.e)**2*(Y.e)**2*C**2 - 2*C*(X.e)**2*(Y.e)**2
    - 2*(X.Y)*(X.e)*(Y.e)*C**2 + 2*(X.Y)*(X.e)*(Y.e)*S**2
W = -2*(X.Y)*(X.e)*(Y.e) + 2*(X.Y)*(X.e)*(Y.e)*C + (X.Y)**2
    + (X.e)**2*(Y.e)**2 + (X.e)**2*(Y.e)**2*C**2 -
     2*C*(X.e)**2*(Y.e)**2
W = {1: -2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2 + (X.e)**2*(Y.e)**2,
     C**2: (X.e)**2*(Y.e)**2,
     C: 2*(X.Y)*(X.e)*(Y.e) - 2*(X.e)**2*(Y.e)**2}

a = (X.e)**2*(Y.e)**2
b = 2*(X.Y)*(X.e)*(Y.e) - 2*(X.e)**2*(Y.e)**2
c = -2*(X.Y)*(X.e)*(Y.e) + (X.Y)**2 + (X.e)**2*(Y.e)**2

Setting to 0 and solving for C gives:
Descriminant D = b^2-4*a*c = 0
C = cosh(alpha) = -b/(2*a) = 1 - (X.Y)/((X.e)*(Y.e))
