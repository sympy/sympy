# Integrands of the form (a+b x+c x**S(2))**p

# Integrands of the form (a+b x+c x**S(2))**p when a=S(0)

# Integrands of the form (b x+c x**S(2))**(p/S(2))

# p>S(0)
[(b*x+c*x**S(2))**(S(7)/S(2)),x,S(6),S(35)/S(6144)*b**S(4)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(3)-S(7)/S(384)*b**S(2)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(5)/S(2))/c**S(2)+S(1)/S(16)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(7)/S(2))/c+S(35)/S(16384)*b**S(8)*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(9)/S(2))-S(35)/S(16384)*b**S(6)*(b+S(2)*c*x)*sqrt(b*x+c*x**S(2))/c**S(4)],

# p<S(0)
[S(1)/(b*x+c*x**S(2))**(S(7)/S(2)),x,S(3),-S(2)/S(5)*(b+S(2)*c*x)/(b**S(2)*(b*x+c*x**S(2))**(S(5)/S(2)))+S(32)/S(15)*c*(b+S(2)*c*x)/(b**S(4)*(b*x+c*x**S(2))**(S(3)/S(2)))-S(256)/S(15)*c**S(2)*(b+S(2)*c*x)/(b**S(6)*sqrt(b*x+c*x**S(2)))],
[S(1)/sqrt(-S(2)*x+x**S(2)),x,S(2),S(2)*arctanh(x/sqrt(-S(2)*x+x**S(2)))],

# Integrands of the form (b x+c x**S(2))**(p/S(3))
[(b*x+c*x**S(2))**(S(4)/S(3)),x,S(6),S(3)/S(55)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(4)/S(3))*(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))/(S(2)**(S(2)/S(3))*c*(-c*(b*x+c*x**S(2))/b**S(2))**(S(4)/S(3)))+S(3)/S(88)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(4)/S(3))*(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(4)/S(3))/(S(2)**(S(2)/S(3))*c*(-c*(b*x+c*x**S(2))/b**S(2))**(S(4)/S(3)))+S(1)/S(55)*S(2)**(S(1)/S(3))*S(3)**(S(3)/S(4))*b**S(2)*(b*x+c*x**S(2))**(S(4)/S(3))*(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3)))*EllipticF((S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))+sqrt(S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3))),sqrt(-S(7)+S(4)*sqrt(S(3))))*sqrt(S(2)-sqrt(S(3)))*sqrt((S(1)+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(2)/S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3)))**S(2))/(c*(b+S(2)*c*x)*(-c*(b*x+c*x**S(2))/b**S(2))**(S(4)/S(3))*sqrt((-S(1)+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3)))**S(2)))],
[(b*x+c*x**S(2))**(S(1)/S(3)),x,S(5),S(3)/S(10)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(1)/S(3))*(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))/(S(2)**(S(2)/S(3))*c*(-c*(b*x+c*x**S(2))/b**S(2))**(S(1)/S(3)))+S(1)/S(5)*S(3)**(S(3)/S(4))*b**S(2)*(b*x+c*x**S(2))**(S(1)/S(3))*(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3)))*EllipticF((S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))+sqrt(S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3))),sqrt(-S(7)+S(4)*sqrt(S(3))))*sqrt(S(2)-sqrt(S(3)))*sqrt((S(1)+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(2)/S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3)))**S(2))/(S(2)**(S(2)/S(3))*c*(b+S(2)*c*x)*(-c*(b*x+c*x**S(2))/b**S(2))**(S(1)/S(3))*sqrt((-S(1)+(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3)))/(S(1)-(S(1)-(b+S(2)*c*x)**S(2)/b**S(2))**(S(1)/S(3))-sqrt(S(3)))**S(2)))],

# Integrands of the form (b x+c x**S(2))**(p/S(4))
[(b*x+c*x**S(2))**(S(5)/S(4)),x,S(5),-S(5)/S(84)*b**S(2)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(1)/S(4))/c**S(2)+S(1)/S(7)*(b+S(2)*c*x)*(b*x+c*x**S(2))**(S(5)/S(4))/c+S(5)/S(84)*b**S(5)*(-c*(b*x+c*x**S(2))/b**S(2))**(S(3)/S(4))*sqrt(cos(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b))**S(2))/cos(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b))*EllipticF(sin(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b)),sqrt(S(2)))/(c**S(3)*(b*x+c*x**S(2))**(S(3)/S(4))*sqrt(S(2)))],
[S(1)/(b*x+c*x**S(2))**(S(13)/S(4)),x,S(6),-S(4)/S(9)*(b+S(2)*c*x)/(b**S(2)*(b*x+c*x**S(2))**(S(9)/S(4)))+S(112)/S(45)*c*(b+S(2)*c*x)/(b**S(4)*(b*x+c*x**S(2))**(S(5)/S(4)))-S(448)/S(15)*c**S(2)*(b+S(2)*c*x)/(b**S(6)*(b*x+c*x**S(2))**(S(1)/S(4)))+S(448)/S(15)*c**S(2)*(-c*(b*x+c*x**S(2))/b**S(2))**(S(1)/S(4))*sqrt(cos(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b))**S(2))/cos(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b))*EllipticE(sin(S(1)/S(2)*arcsin(S(1)+S(2)*c*x/b)),sqrt(S(2)))*sqrt(S(2))/(b**S(5)*(b*x+c*x**S(2))**(S(1)/S(4)))],

# Integrands of the form (b x+c x**S(2))**p with p symbolic
[(b*x+c*x**S(2))**p,x,S(1),-(-c*x/b)**(-S(1)-p)*(b*x+c*x**S(2))**(S(1)+p)*hypergeom([-p,S(1)+p],[S(2)+p],(b+c*x)/b)/(b*(S(1)+p))],

# Integrands of the form (a+b x+c x**S(2))**p when b**S(2)-S(4) a c=S(0)

# Integrands of the form (a**S(2)+S(2) a b x+b**S(2) x**S(2))**(p/S(2))

[sqrt(S(4)+S(12)*x+S(9)*x**S(2)),x,S(1),S(1)/S(6)*(S(2)+S(3)*x)*sqrt(S(4)+S(12)*x+S(9)*x**S(2))],
[S(1)/sqrt(-S(4)+S(12)*x-S(9)*x**S(2)),x,S(2),-S(1)/S(3)*(S(2)-S(3)*x)*log(S(2)-S(3)*x)/sqrt(-S(4)+S(12)*x-S(9)*x**S(2))],
[S(1)/sqrt(-S(4)-S(12)*x-S(9)*x**S(2)),x,S(2),S(1)/S(3)*(S(2)+S(3)*x)*log(S(2)+S(3)*x)/sqrt(-S(4)-S(12)*x-S(9)*x**S(2))],

# Integrands of the form (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p with p symbolic
[(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p,x,S(1),(a+b*x)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b*(S(1)+S(2)*p))],

# Integrands of the form (a+b x+c x**S(2))**p

# Integrands of the form (a+b x+c x**S(2))**p

[(S(1)/S(4)*(-S(1)+b**S(2))/c+b*x+c*x**S(2))**S(5),x,S(3),S(1)/S(384)*(S(1)-b-S(2)*c*x)**S(6)/c**S(6)-S(5)/S(896)*(S(1)-b-S(2)*c*x)**S(7)/c**S(6)+S(5)/S(1024)*(S(1)-b-S(2)*c*x)**S(8)/c**S(6)-S(5)/S(2304)*(S(1)-b-S(2)*c*x)**S(9)/c**S(6)+S(1)/S(2048)*(S(1)-b-S(2)*c*x)**S(10)/c**S(6)-S(1)/S(22528)*(S(1)-b-S(2)*c*x)**S(11)/c**S(6)],
[S(1)/(S(1)+x**S(2)+S(2)*x*cos(S(1)/S(7)*Pi)),x,S(2),arctan((x+cos(S(1)/S(7)*Pi))*csc(S(1)/S(7)*Pi))*csc(S(1)/S(7)*Pi)],

# Integrands of the form (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[sqrt(S(5)-S(6)*x+S(9)*x**S(2)),x,S(3),S(2)/S(3)*arcsinh(S(1)/S(2)*(-S(1)+S(3)*x))-S(1)/S(6)*(S(1)-S(3)*x)*sqrt(S(5)-S(6)*x+S(9)*x**S(2))],

# p<S(0)
[S(1)/sqrt(S(5)-S(6)*x+S(9)*x**S(2)),x,S(2),S(1)/S(3)*arcsinh(S(1)/S(2)*(-S(1)+S(3)*x))],
[S(1)/sqrt(S(3)-S(4)*x-S(4)*x**S(2)),x,S(2),S(1)/S(2)*arcsin(S(1)/S(2)+x),-S(1)/S(2)*arcsin(S(1)/S(2)*(-S(1)-S(2)*x))],

# Integrands of the form (a+b x+c x**S(2))**p with p symbolic
[(a+b*x+c*x**S(2))**p,x,S(1),-S(2)**(S(1)+p)*(a+b*x+c*x**S(2))**(S(1)+p)*hypergeom([-p,S(1)+p],[S(2)+p],S(1)/S(2)*(b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))*((-b-S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))**(-S(1)-p)/((S(1)+p)*sqrt(b**S(2)-S(4)*a*c))],
[(S(3)+S(4)*x-S(5)*x**S(2))**p,x,S(2),-S(5)**(-S(1)-p)*S(19)**p*(S(2)-S(5)*x)*hypergeom([S(1)/S(2),-p],[S(3)/S(2)],S(1)/S(19)*(S(2)-S(5)*x)**S(2))],

# Integrands of the form (d x)**m (a+b x+c x**S(2))**p

# Integrands of the form (d x)**m (a+b x+c x**S(2))**p when a=S(0)

# Integrands of the form x**m (b x+c x**S(2))**p

[x**S(2)*(b*x+c*x**S(2)),x,S(2),S(1)/S(4)*b*x**S(4)+S(1)/S(5)*c*x**S(5)],
[S(1)/(x*(b*x+c*x**S(2))**S(3)),x,S(3),(-S(1)/S(3))/(b**S(3)*x**S(3))+S(3)/S(2)*c/(b**S(4)*x**S(2))-S(6)*c**S(2)/(b**S(5)*x)-S(1)/S(2)*c**S(3)/(b**S(4)*(b+c*x)**S(2))-S(4)*c**S(3)/(b**S(5)*(b+c*x))-S(10)*c**S(3)*log(x)/b**S(6)+S(10)*c**S(3)*log(b+c*x)/b**S(6)],

# Integrands of the form x**m (b x+c x**S(2))**(p/S(2))

# p>S(0)
[x**S(3)*(b*x+c*x**S(2))**(S(1)/S(2)),x,S(6),S(7)/S(48)*b**S(2)*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(3)-S(7)/S(40)*b*x*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(5)*x**S(2)*(b*x+c*x**S(2))**(S(3)/S(2))/c+S(7)/S(128)*b**S(5)*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(9)/S(2))-S(7)/S(128)*b**S(3)*(b+S(2)*c*x)*sqrt(b*x+c*x**S(2))/c**S(4)],
[x**S(2)*(b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),-S(5)/S(24)*b*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(4)*x*(b*x+c*x**S(2))**(S(3)/S(2))/c-S(5)/S(64)*b**S(4)*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(7)/S(2))+S(5)/S(64)*b**S(2)*(b+S(2)*c*x)*sqrt(b*x+c*x**S(2))/c**S(3)],

# p<S(0)
[x**S(4)/(b*x+c*x**S(2))**(S(1)/S(2)),x,S(6),S(35)/S(64)*b**S(4)*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(9)/S(2))-S(35)/S(64)*b**S(3)*sqrt(b*x+c*x**S(2))/c**S(4)+S(35)/S(96)*b**S(2)*x*sqrt(b*x+c*x**S(2))/c**S(3)-S(7)/S(24)*b*x**S(2)*sqrt(b*x+c*x**S(2))/c**S(2)+S(1)/S(4)*x**S(3)*sqrt(b*x+c*x**S(2))/c],
[x**S(3)/(b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),-S(5)/S(8)*b**S(3)*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(7)/S(2))+S(5)/S(8)*b**S(2)*sqrt(b*x+c*x**S(2))/c**S(3)-S(5)/S(12)*b*x*sqrt(b*x+c*x**S(2))/c**S(2)+S(1)/S(3)*x**S(2)*sqrt(b*x+c*x**S(2))/c],

# Integrands of the form x**(m/S(2)) (b x+c x**S(2))**p

# p>S(0)
[x**(S(7)/S(2))*(b*x+c*x**S(2)),x,S(2),S(2)/S(11)*b*x**(S(11)/S(2))+S(2)/S(13)*c*x**(S(13)/S(2))],
[x**(S(5)/S(2))*(b*x+c*x**S(2)),x,S(2),S(2)/S(9)*b*x**(S(9)/S(2))+S(2)/S(11)*c*x**(S(11)/S(2))],

# p<S(0)
[x**(S(7)/S(2))/(b*x+c*x**S(2)),x,S(6),-S(2)/S(3)*b*x**(S(3)/S(2))/c**S(2)+S(2)/S(5)*x**(S(5)/S(2))/c-S(2)*b**(S(5)/S(2))*arctan(sqrt(c)*sqrt(x)/sqrt(b))/c**(S(7)/S(2))+S(2)*b**S(2)*sqrt(x)/c**S(3)],
[x**(S(5)/S(2))/(b*x+c*x**S(2)),x,S(5),S(2)/S(3)*x**(S(3)/S(2))/c+S(2)*b**(S(3)/S(2))*arctan(sqrt(c)*sqrt(x)/sqrt(b))/c**(S(5)/S(2))-S(2)*b*sqrt(x)/c**S(2)],

# Integrands of the form x**(m/S(2)) (b x+c x**S(2))**(p/S(2))

# p>S(0)
[x**(S(7)/S(2))*(b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),S(256)/S(3465)*b**S(4)*(b*x+c*x**S(2))**(S(3)/S(2))/(c**S(5)*x**(S(3)/S(2)))-S(16)/S(99)*b*x**(S(3)/S(2))*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(2)/S(11)*x**(S(5)/S(2))*(b*x+c*x**S(2))**(S(3)/S(2))/c-S(128)/S(1155)*b**S(3)*(b*x+c*x**S(2))**(S(3)/S(2))/(c**S(4)*sqrt(x))+S(32)/S(231)*b**S(2)*(b*x+c*x**S(2))**(S(3)/S(2))*sqrt(x)/c**S(3)],
[S(1)/(x**(S(3)/S(2))*(b*x+c*x**S(2))**(S(3)/S(2))),x,S(5),-S(15)/S(4)*c**S(2)*arctanh(sqrt(b*x+c*x**S(2))/(sqrt(b)*sqrt(x)))/b**(S(7)/S(2))+(-S(1)/S(2))/(b*x**(S(3)/S(2))*sqrt(b*x+c*x**S(2)))+S(5)/S(4)*c/(b**S(2)*sqrt(x)*sqrt(b*x+c*x**S(2)))+S(15)/S(4)*c**S(2)*sqrt(x)/(b**S(3)*sqrt(b*x+c*x**S(2)))],
[S(1)/(x**(S(5)/S(2))*(b*x+c*x**S(2))**(S(3)/S(2))),x,S(6),S(35)/S(8)*c**S(3)*arctanh(sqrt(b*x+c*x**S(2))/(sqrt(b)*sqrt(x)))/b**(S(9)/S(2))+(-S(1)/S(3))/(b*x**(S(5)/S(2))*sqrt(b*x+c*x**S(2)))+S(7)/S(12)*c/(b**S(2)*x**(S(3)/S(2))*sqrt(b*x+c*x**S(2)))-S(35)/S(24)*c**S(2)/(b**S(3)*sqrt(x)*sqrt(b*x+c*x**S(2)))-S(35)/S(8)*c**S(3)*sqrt(x)/(b**S(4)*sqrt(b*x+c*x**S(2)))],

# Integrands of the form (d x)**m (b x+c x**S(2))**p with m symbolic
[(d*x)**m*(b*x+c*x**S(2))**S(3),x,S(4),b**S(3)*(d*x)**(S(4)+m)/(d**S(4)*(S(4)+m))+S(3)*b**S(2)*c*(d*x)**(S(5)+m)/(d**S(5)*(S(5)+m))+S(3)*b*c**S(2)*(d*x)**(S(6)+m)/(d**S(6)*(S(6)+m))+c**S(3)*(d*x)**(S(7)+m)/(d**S(7)*(S(7)+m))],
[(d*x)**m*(b*x+c*x**S(2))**S(2),x,S(4),b**S(2)*(d*x)**(S(3)+m)/(d**S(3)*(S(3)+m))+S(2)*b*c*(d*x)**(S(4)+m)/(d**S(4)*(S(4)+m))+c**S(2)*(d*x)**(S(5)+m)/(d**S(5)*(S(5)+m))],

# Integrands of the form (d x)**m (b x+c x**S(2))**p with p symbolic
[(d*x)**m*(b*x+c*x**S(2))**p,x,S(3),x*(d*x)**m*(b*x+c*x**S(2))**p*hypergeom([-p,S(1)+m+p],[S(2)+m+p],-c*x/b)/((S(1)+m+p)*(S(1)+c*x/b)**p)],
[x**S(3)*(b*x+c*x**S(2))**p,x,S(3),x**S(4)*(b*x+c*x**S(2))**p*hypergeom([-p,S(4)+p],[S(5)+p],-c*x/b)/((S(4)+p)*(S(1)+c*x/b)**p)],

# Integrands of the form (d x)**m (a+b x+c x**S(2))**p when b**S(2)-S(4) a c=S(0)

# Integrands of the form x**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p

# p>S(0)
[x**S(4)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(1)/S(5)*a**S(2)*x**S(5)+S(1)/S(3)*a*b*x**S(6)+S(1)/S(7)*b**S(2)*x**S(7)],
[x**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(1)/S(4)*a**S(2)*x**S(4)+S(2)/S(5)*a*b*x**S(5)+S(1)/S(6)*b**S(2)*x**S(6)],

#  {(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2))**S(2)/x**S(7), x, S(3), -(a + b*x)**S(5)/(S(6)*a*x**S(6)) + (b*(a + b*x)**S(5))/(S(30)*a**S(2)*x**S(5))} 
[(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(2)/x**S(7),x,S(3),-S(1)/S(6)*a**S(4)/x**S(6)-S(4)/S(5)*a**S(3)*b/x**S(5)-S(3)/S(2)*a**S(2)*b**S(2)/x**S(4)-S(4)/S(3)*a*b**S(3)/x**S(3)-S(1)/S(2)*b**S(4)/x**S(2)],
[S(1)/(x**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(3)),x,S(3),(-S(1))/(a**S(6)*x)-S(1)/S(5)*b/(a**S(2)*(a+b*x)**S(5))-S(1)/S(2)*b/(a**S(3)*(a+b*x)**S(4))-b/(a**S(4)*(a+b*x)**S(3))-S(2)*b/(a**S(5)*(a+b*x)**S(2))-S(5)*b/(a**S(6)*(a+b*x))-S(6)*b*log(x)/a**S(7)+S(6)*b*log(a+b*x)/a**S(7)],
[S(1)/(x**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(3)),x,S(3),(-S(1)/S(2))/(a**S(6)*x**S(2))+S(6)*b/(a**S(7)*x)+S(1)/S(5)*b**S(2)/(a**S(3)*(a+b*x)**S(5))+S(3)/S(4)*b**S(2)/(a**S(4)*(a+b*x)**S(4))+S(2)*b**S(2)/(a**S(5)*(a+b*x)**S(3))+S(5)*b**S(2)/(a**S(6)*(a+b*x)**S(2))+S(15)*b**S(2)/(a**S(7)*(a+b*x))+S(21)*b**S(2)*log(x)/a**S(8)-S(21)*b**S(2)*log(a+b*x)/a**S(8)],

# Integrands of the form x**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**(p/S(2))

[x**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(1)/S(6)*x**S(5)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))+S(1)/S(30)*a*x**S(5)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(a+b*x)],
[x**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(1)/S(5)*x**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))+S(1)/S(20)*a*x**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(a+b*x)],
[x/sqrt(-S(4)-S(12)*x-S(9)*x**S(2)),x,S(3),-S(2)/S(9)*(S(2)+S(3)*x)*log(S(2)+S(3)*x)/sqrt(-S(4)-S(12)*x-S(9)*x**S(2))-S(1)/S(9)*sqrt(-S(4)-S(12)*x-S(9)*x**S(2))],

# Integrands of the form x**(m/S(2)) (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p

# p>S(0)
[x**(S(7)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(2)/S(9)*a**S(2)*x**(S(9)/S(2))+S(4)/S(11)*a*b*x**(S(11)/S(2))+S(2)/S(13)*b**S(2)*x**(S(13)/S(2))],
[x**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(2)/S(7)*a**S(2)*x**(S(7)/S(2))+S(4)/S(9)*a*b*x**(S(9)/S(2))+S(2)/S(11)*b**S(2)*x**(S(11)/S(2))],

# p<S(0)
[x**(S(9)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(8),S(3)*a**S(2)*x**(S(3)/S(2))/b**S(4)-S(9)/S(5)*a*x**(S(5)/S(2))/b**S(3)+S(9)/S(7)*x**(S(7)/S(2))/b**S(2)-x**(S(9)/S(2))/(b*(a+b*x))+S(9)*a**(S(7)/S(2))*arctan(sqrt(b)*sqrt(x)/sqrt(a))/b**(S(11)/S(2))-S(9)*a**S(3)*sqrt(x)/b**S(5)],
[S(1)/(x**(S(7)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(3)),x,S(11),(-S(9009)/S(640))/(a**S(6)*x**(S(5)/S(2)))+S(3003)/S(128)*b/(a**S(7)*x**(S(3)/S(2)))+S(1)/S(5)/(a*x**(S(5)/S(2))*(a+b*x)**S(5))+S(3)/S(8)/(a**S(2)*x**(S(5)/S(2))*(a+b*x)**S(4))+S(13)/S(16)/(a**S(3)*x**(S(5)/S(2))*(a+b*x)**S(3))+S(143)/S(64)/(a**S(4)*x**(S(5)/S(2))*(a+b*x)**S(2))+S(1287)/S(128)/(a**S(5)*x**(S(5)/S(2))*(a+b*x))-S(9009)/S(128)*b**(S(5)/S(2))*arctan(sqrt(b)*sqrt(x)/sqrt(a))/a**(S(17)/S(2))-S(9009)/S(128)*b**S(2)/(a**S(8)*sqrt(x))],

# Integrands of the form x**(m/S(2)) (a**S(2)+S(2) a b x+b**S(2) x**S(2))**(p/S(2))

# p>S(0)
[x**(S(7)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(2),S(2)/S(11)*x**(S(9)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))+S(4)/S(99)*a*x**(S(9)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(a+b*x)],
[(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(5)/S(2))/x**(S(7)/S(2)),x,S(4),-S(4)/S(3)*b*(a+b*x)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(3)/S(2))/x**(S(3)/S(2))-S(2)/S(5)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(5)/S(2))/x**(S(5)/S(2))-S(32)/S(3)*b**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(3)/S(2))/sqrt(x)+S(256)/S(15)*a*b**S(3)*sqrt(x)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))+S(512)/S(15)*a**S(2)*b**S(3)*sqrt(x)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(a+b*x)+S(64)/S(5)*b**S(3)*(a+b*x)*sqrt(x)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))],

# p<S(0)
[x**(S(7)/S(2))/sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(7),S(2)/S(3)*a**S(2)*x**(S(3)/S(2))*(a+b*x)/(b**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(2)/S(5)*a*x**(S(5)/S(2))*(a+b*x)/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)/S(7)*x**(S(7)/S(2))*(a+b*x)/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)*a**(S(7)/S(2))*(a+b*x)*arctan(sqrt(b)*sqrt(x)/sqrt(a))/(b**(S(9)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(2)*a**S(3)*(a+b*x)*sqrt(x)/(b**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],
[S(1)/(x**(S(7)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(5)/S(2))),x,S(8),S(1)/S(4)*(a+b*x)/(a*x**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(5)/S(2)))+S(13)/S(24)/(a**S(2)*x**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(3)/S(2)))+S(143)/S(96)*(a+b*x)/(a**S(3)*x**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(3)/S(2)))+S(429)/S(64)/(a**S(4)*x**(S(5)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(3003)/S(320)*(a+b*x)/(a**S(5)*x**(S(5)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(1001)/S(64)*b*(a+b*x)/(a**S(6)*x**(S(3)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(3003)/S(64)*b**(S(5)/S(2))*(a+b*x)*arctan(sqrt(b)*sqrt(x)/sqrt(a))/(a**(S(15)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(3003)/S(64)*b**S(2)*(a+b*x)/(a**S(7)*sqrt(x)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],

# Integrands of the form (d x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p with m symbolic
[(d*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(3),x,S(3),a**S(6)*(d*x)**(S(1)+m)/(d*(S(1)+m))+S(6)*a**S(5)*b*(d*x)**(S(2)+m)/(d**S(2)*(S(2)+m))+S(15)*a**S(4)*b**S(2)*(d*x)**(S(3)+m)/(d**S(3)*(S(3)+m))+S(20)*a**S(3)*b**S(3)*(d*x)**(S(4)+m)/(d**S(4)*(S(4)+m))+S(15)*a**S(2)*b**S(4)*(d*x)**(S(5)+m)/(d**S(5)*(S(5)+m))+S(6)*a*b**S(5)*(d*x)**(S(6)+m)/(d**S(6)*(S(6)+m))+b**S(6)*(d*x)**(S(7)+m)/(d**S(7)*(S(7)+m))],
[(d*x)**m/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(5)/S(2)),x,S(2),(d*x)**(S(1)+m)*(a+b*x)*hypergeom([S(5),S(1)+m],[S(2)+m],-b*x/a)/(a**S(5)*d*(S(1)+m)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],

# Integrands of the form (d x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p with p symbolic
[(d*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p,x,S(3),(d*x)**(S(1)+m)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p*hypergeom([S(1)+m,-S(2)*p],[S(2)+m],-b*x/a)/(d*(S(1)+m)*(S(1)+b*x/a)**(S(2)*p))],

# Integrands of the form (d x)**m (a+b x+c x**S(2))**p

# Integrands of the form x**m (a+b x+c x**S(2))**p

# p>S(0)

#  These should ALL require only S(2) rule applications! 
[x**m*(S(3)-S(4)*x+x**S(2))**S(2),x,S(2),S(9)*x**(S(1)+m)/(S(1)+m)-S(24)*x**(S(2)+m)/(S(2)+m)+S(22)*x**(S(3)+m)/(S(3)+m)-S(8)*x**(S(4)+m)/(S(4)+m)+x**(S(5)+m)/(S(5)+m)],

# p<S(0)
[x**S(7)/(a+b*x+c*x**S(2))**S(3),x,S(9),-S(3)*b*(S(2)*b**S(2)-S(9)*a*c)*(b**S(2)-S(3)*a*c)*x/(c**S(4)*(b**S(2)-S(4)*a*c)**S(2))+S(3)/S(2)*(S(2)*b**S(4)-S(13)*a*b**S(2)*c+S(16)*a**S(2)*c**S(2))*x**S(2)/(c**S(3)*(b**S(2)-S(4)*a*c)**S(2))-b*(S(2)*b**S(2)-S(11)*a*c)*x**S(3)/(c**S(2)*(b**S(2)-S(4)*a*c)**S(2))+S(1)/S(2)*x**S(6)*(S(2)*a+b*x)/((b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**S(2))+S(3)/S(2)*x**S(4)*(a*(b**S(2)-S(8)*a*c)+b*(b**S(2)-S(6)*a*c)*x)/(c*(b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2)))+S(3)*b*(S(2)*b**S(6)-S(21)*a*b**S(4)*c+S(70)*a**S(2)*b**S(2)*c**S(2)-S(70)*a**S(3)*c**S(3))*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c**S(5)*(b**S(2)-S(4)*a*c)**(S(5)/S(2)))+S(3)/S(2)*(S(2)*b**S(2)-a*c)*log(a+b*x+c*x**S(2))/c**S(5)],
[x**S(6)/(a+b*x+c*x**S(2))**S(3),x,S(9),S(3)*(b**S(4)-S(7)*a*b**S(2)*c+S(10)*a**S(2)*c**S(2))*x/(c**S(3)*(b**S(2)-S(4)*a*c)**S(2))-S(3)/S(2)*b*(b**S(2)-S(6)*a*c)*x**S(2)/(c**S(2)*(b**S(2)-S(4)*a*c)**S(2))+S(1)/S(2)*x**S(5)*(S(2)*a+b*x)/((b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**S(2))+x**S(3)*(a*(b**S(2)-S(10)*a*c)+b*(b**S(2)-S(7)*a*c)*x)/(c*(b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2)))-S(3)*(b**S(6)-S(10)*a*b**S(4)*c+S(30)*a**S(2)*b**S(2)*c**S(2)-S(20)*a**S(3)*c**S(3))*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c**S(4)*(b**S(2)-S(4)*a*c)**(S(5)/S(2)))-S(3)/S(2)*b*log(a+b*x+c*x**S(2))/c**S(4)],

# Integrands of the form (d x)**(m/S(2)) (a+b x+c x**S(2))**p

[(a+b*x+c*x**S(2))*sqrt(d*x),x,S(2),S(2)/S(3)*a*(d*x)**(S(3)/S(2))/d+S(2)/S(5)*b*(d*x)**(S(5)/S(2))/d**S(2)+S(2)/S(7)*c*(d*x)**(S(7)/S(2))/d**S(3)],
[x**(S(9)/S(2))/(a+b*x+c*x**S(2))**S(3),x,S(7),S(1)/S(2)*x**(S(7)/S(2))*(S(2)*a+b*x)/((b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**S(2))+S(1)/S(4)*x**(S(3)/S(2))*(a*(b**S(2)-S(28)*a*c)+b*(b**S(2)-S(16)*a*c)*x)/(c*(b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2)))-S(3)/S(4)*b*(b**S(2)-S(8)*a*c)*sqrt(x)/(c**S(2)*(b**S(2)-S(4)*a*c)**S(2))+S(3)/S(4)*arctan(sqrt(S(2))*sqrt(c)*sqrt(x)/sqrt(b-sqrt(b**S(2)-S(4)*a*c)))*(b**S(4)-S(9)*a*b**S(2)*c+S(28)*a**S(2)*c**S(2)+(-b**S(5)+S(11)*a*b**S(3)*c-S(44)*a**S(2)*b*c**S(2))/sqrt(b**S(2)-S(4)*a*c))/(c**(S(5)/S(2))*(b**S(2)-S(4)*a*c)**S(2)*sqrt(S(2))*sqrt(b-sqrt(b**S(2)-S(4)*a*c)))+S(3)/S(4)*arctan(sqrt(S(2))*sqrt(c)*sqrt(x)/sqrt(b+sqrt(b**S(2)-S(4)*a*c)))*(b**S(4)-S(9)*a*b**S(2)*c+S(28)*a**S(2)*c**S(2)+(b**S(5)-S(11)*a*b**S(3)*c+S(44)*a**S(2)*b*c**S(2))/sqrt(b**S(2)-S(4)*a*c))/(c**(S(5)/S(2))*(b**S(2)-S(4)*a*c)**S(2)*sqrt(S(2))*sqrt(b+sqrt(b**S(2)-S(4)*a*c)))],
[S(1)/(x**(S(3)/S(2))*(a+b*x+c*x**S(2))**S(3)),x,S(7),-S(3)/S(4)*(S(5)*b**S(2)-S(12)*a*c)*(b**S(2)-S(5)*a*c)/(a**S(3)*(b**S(2)-S(4)*a*c)**S(2)*sqrt(x))+S(1)/S(2)*(b**S(2)-S(2)*a*c+b*c*x)/(a*(b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**S(2)*sqrt(x))+S(1)/S(4)*(S(5)*b**S(4)-S(35)*a*b**S(2)*c+S(36)*a**S(2)*c**S(2)+b*c*(S(5)*b**S(2)-S(32)*a*c)*x)/(a**S(2)*(b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2))*sqrt(x))-S(3)/S(4)*arctan(sqrt(S(2))*sqrt(c)*sqrt(x)/sqrt(b-sqrt(b**S(2)-S(4)*a*c)))*sqrt(c)*(S(5)*b**S(5)-S(47)*a*b**S(3)*c+S(124)*a**S(2)*b*c**S(2)+(S(5)*b**S(4)-S(37)*a*b**S(2)*c+S(60)*a**S(2)*c**S(2))*sqrt(b**S(2)-S(4)*a*c))/(a**S(3)*(b**S(2)-S(4)*a*c)**(S(5)/S(2))*sqrt(S(2))*sqrt(b-sqrt(b**S(2)-S(4)*a*c)))+S(3)/S(4)*arctan(sqrt(S(2))*sqrt(c)*sqrt(x)/sqrt(b+sqrt(b**S(2)-S(4)*a*c)))*sqrt(c)*(S(5)*b**S(5)-S(47)*a*b**S(3)*c+S(124)*a**S(2)*b*c**S(2)-S(5)*b**S(4)*sqrt(b**S(2)-S(4)*a*c)+S(37)*a*b**S(2)*c*sqrt(b**S(2)-S(4)*a*c)-S(60)*a**S(2)*c**S(2)*sqrt(b**S(2)-S(4)*a*c))/(a**S(3)*(b**S(2)-S(4)*a*c)**(S(5)/S(2))*sqrt(S(2))*sqrt(b+sqrt(b**S(2)-S(4)*a*c)))],

# Integrands of the form x**m (a+b x+c x**S(2))**(p/S(2))
[S(1)/(x**S(3)*sqrt(S(1)+x+x**S(2))),x,S(4),S(1)/S(8)*arctanh(S(1)/S(2)*(S(2)+x)/sqrt(S(1)+x+x**S(2)))-S(1)/S(2)*sqrt(S(1)+x+x**S(2))/x**S(2)+S(3)/S(4)*sqrt(S(1)+x+x**S(2))/x],
[S(1)/x+(-S(1))/(x*sqrt(S(1)+b*x+c*x**S(2))),x,S(3),log(-S(2)-b*x-S(2)*sqrt(S(1)+b*x+c*x**S(2))),arctanh(S(1)/S(2)*(S(2)+b*x)/sqrt(S(1)+b*x+c*x**S(2)))+log(x)],

# Integrands of the form (d x)**(m/S(2)) (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[(d*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2)),x,S(8),S(2)/S(9)*d*(d*x)**(S(3)/S(2))*(a+b*x+c*x**S(2))**(S(3)/S(2))/c-S(4)/S(21)*b*d**S(2)*(a+b*x+c*x**S(2))**(S(3)/S(2))*sqrt(d*x)/c**S(2)-S(4)/S(315)*(S(8)*b**S(4)-S(36)*a*b**S(2)*c+S(21)*a**S(2)*c**S(2))*d**S(3)*x*sqrt(a+b*x+c*x**S(2))/(c**(S(7)/S(2))*(sqrt(a)+x*sqrt(c))*sqrt(d*x))+S(2)/S(315)*d**S(2)*(b*(S(8)*b**S(2)+S(3)*a*c)+S(3)*c*(S(8)*b**S(2)-S(7)*a*c)*x)*sqrt(d*x)*sqrt(a+b*x+c*x**S(2))/c**S(3)+S(4)/S(315)*a**(S(1)/S(4))*(S(8)*b**S(4)-S(36)*a*b**S(2)*c+S(21)*a**S(2)*c**S(2))*d**S(3)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(15)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))-S(1)/S(315)*a**(S(1)/S(4))*d**S(3)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*(S(16)*b**S(4)-S(72)*a*b**S(2)*c+S(42)*a**S(2)*c**S(2)+b*(S(8)*b**S(2)-S(27)*a*c)*sqrt(a)*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(15)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))],

# p<S(0)
[(d*x)**(S(7)/S(2))/sqrt(a+b*x+c*x**S(2)),x,S(8),-S(12)/S(35)*b*d**S(2)*(d*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2))/c**S(2)+S(2)/S(7)*d*(d*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2))/c-S(16)/S(105)*b*(S(6)*b**S(2)-S(13)*a*c)*d**S(4)*x*sqrt(a+b*x+c*x**S(2))/(c**(S(7)/S(2))*(sqrt(a)+x*sqrt(c))*sqrt(d*x))+S(2)/S(105)*(S(24)*b**S(2)-S(25)*a*c)*d**S(3)*sqrt(d*x)*sqrt(a+b*x+c*x**S(2))/c**S(3)+S(16)/S(105)*a**(S(1)/S(4))*b*(S(6)*b**S(2)-S(13)*a*c)*d**S(4)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(15)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))-S(1)/S(105)*a**(S(1)/S(4))*d**S(4)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*(S(8)*b*(S(6)*b**S(2)-S(13)*a*c)+(S(24)*b**S(2)-S(25)*a*c)*sqrt(a)*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(15)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))],
[(d*x)**(S(5)/S(2))/sqrt(a+b*x+c*x**S(2)),x,S(7),S(2)/S(5)*d*(d*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2))/c+S(2)/S(15)*(S(8)*b**S(2)-S(9)*a*c)*d**S(3)*x*sqrt(a+b*x+c*x**S(2))/(c**(S(5)/S(2))*(sqrt(a)+x*sqrt(c))*sqrt(d*x))-S(8)/S(15)*b*d**S(2)*sqrt(d*x)*sqrt(a+b*x+c*x**S(2))/c**S(2)-S(2)/S(15)*a**(S(1)/S(4))*(S(8)*b**S(2)-S(9)*a*c)*d**S(3)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(11)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))+S(1)/S(15)*a**(S(1)/S(4))*d**S(3)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*sqrt(x)/a**(S(1)/S(4)))),sqrt(S(1)/S(4)*(S(2)-b/(sqrt(a)*sqrt(c)))))*(sqrt(a)+x*sqrt(c))*(S(8)*b**S(2)-S(9)*a*c+S(4)*b*sqrt(a)*sqrt(c))*sqrt(x)*sqrt((a+b*x+c*x**S(2))/(sqrt(a)+x*sqrt(c))**S(2))/(c**(S(11)/S(4))*sqrt(d*x)*sqrt(a+b*x+c*x**S(2)))],

# Integrands of the form (d x)**m (a+b x+c x**S(2))**p with p symbolic
[(d*x)**m*(a+b*x+c*x**S(2))**p,x,S(2),(d*x)**(S(1)+m)*(a+b*x+c*x**S(2))**p*AppellFS(1)(S(1)+m,-p,-p,S(2)+m,-S(2)*c*x/(b-sqrt(b**S(2)-S(4)*a*c)),-S(2)*c*x/(b+sqrt(b**S(2)-S(4)*a*c)))/(d*(S(1)+m)*(S(1)+S(2)*c*x/(b-sqrt(b**S(2)-S(4)*a*c)))**p*(S(1)+S(2)*c*x/(b+sqrt(b**S(2)-S(4)*a*c)))**p)],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when a=S(0)

# Integrands of the form (b d+S(2) c d x)**m (b x+c x**S(2))**p

# p>S(0)

# p<S(0)
[(S(1)+x)/(S(2)*x+x**S(2)),x,S(2),S(1)/S(2)*log(S(2)*x+x**S(2))],
[(a+S(2)*b*x)/(a*x+b*x**S(2)),x,S(2),log(a*x+b*x**S(2))],

# Integrands of the form (d+e x)**m (b x+c x**S(2))**p

# p>S(0)
[(d+e*x)**S(4)*(b*x+c*x**S(2)),x,S(3),S(1)/S(5)*d*(c*d-b*e)*(d+e*x)**S(5)/e**S(3)-S(1)/S(6)*(S(2)*c*d-b*e)*(d+e*x)**S(6)/e**S(3)+S(1)/S(7)*c*(d+e*x)**S(7)/e**S(3)],
[(d+e*x)**S(3)*(b*x+c*x**S(2)),x,S(3),S(1)/S(4)*d*(c*d-b*e)*(d+e*x)**S(4)/e**S(3)-S(1)/S(5)*(S(2)*c*d-b*e)*(d+e*x)**S(5)/e**S(3)+S(1)/S(6)*c*(d+e*x)**S(6)/e**S(3)],

# p<S(0)
[(d+e*x)**S(4)/(b*x+c*x**S(2)),x,S(3),e**S(2)*(S(6)*c**S(2)*d**S(2)-S(4)*b*c*d*e+b**S(2)*e**S(2))*x/c**S(3)+S(1)/S(2)*e**S(3)*(S(4)*c*d-b*e)*x**S(2)/c**S(2)+S(1)/S(3)*e**S(4)*x**S(3)/c+d**S(4)*log(x)/b-(c*d-b*e)**S(4)*log(b+c*x)/(b*c**S(4))],
[(d+e*x)**S(3)/(b*x+c*x**S(2)),x,S(3),e**S(2)*(S(3)*c*d-b*e)*x/c**S(2)+S(1)/S(2)*e**S(3)*x**S(2)/c+d**S(3)*log(x)/b-(c*d-b*e)**S(3)*log(b+c*x)/(b*c**S(3))],

# Integrands of the form (d+e x)**m (b x+c x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**S(3)*(b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),S(1)/S(5)*e*(d+e*x)**S(2)*(b*x+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(240)*e*(S(192)*c**S(2)*d**S(2)-S(150)*b*c*d*e+S(35)*b**S(2)*e**S(2)+S(42)*c*e*(S(2)*c*d-b*e)*x)*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(3)-S(1)/S(128)*b**S(2)*(S(2)*c*d-b*e)*(S(16)*c**S(2)*d**S(2)-S(16)*b*c*d*e+S(7)*b**S(2)*e**S(2))*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(9)/S(2))+S(1)/S(128)*(S(2)*c*d-b*e)*(S(16)*c**S(2)*d**S(2)-S(16)*b*c*d*e+S(7)*b**S(2)*e**S(2))*(b+S(2)*c*x)*sqrt(b*x+c*x**S(2))/c**S(4)],
[(d+e*x)**S(2)*(b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),S(5)/S(24)*e*(S(2)*c*d-b*e)*(b*x+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(4)*e*(d+e*x)*(b*x+c*x**S(2))**(S(3)/S(2))/c-S(1)/S(64)*b**S(2)*(S(16)*c**S(2)*d**S(2)-S(16)*b*c*d*e+S(5)*b**S(2)*e**S(2))*arctanh(x*sqrt(c)/sqrt(b*x+c*x**S(2)))/c**(S(7)/S(2))+S(1)/S(64)*(S(16)*c**S(2)*d**S(2)-S(16)*b*c*d*e+S(5)*b**S(2)*e**S(2))*(b+S(2)*c*x)*sqrt(b*x+c*x**S(2))/c**S(3)],

#  Integrands of the form (d+e*x)**m*(b*x+c*x**S(2))**p where S(2)*c*d-b*e=S(0) 
[sqrt(S(2)*x+x**S(2))/(S(1)+x),x,S(3),-arctan(sqrt(S(2)*x+x**S(2)))+sqrt(S(2)*x+x**S(2))],
[S(1)/((S(2)+x)*sqrt(S(2)*x+x**S(2))),x,S(1),sqrt(S(2)*x+x**S(2))/(S(2)+x)],

# Integrands of the form (d+e x)**(m/S(2)) (b x+c x**S(2))**p

[(d+e*x)**(S(7)/S(2))*(b*x+c*x**S(2)),x,S(3),S(2)/S(9)*d*(c*d-b*e)*(d+e*x)**(S(9)/S(2))/e**S(3)-S(2)/S(11)*(S(2)*c*d-b*e)*(d+e*x)**(S(11)/S(2))/e**S(3)+S(2)/S(13)*c*(d+e*x)**(S(13)/S(2))/e**S(3)],
[(d+e*x)**(S(5)/S(2))*(b*x+c*x**S(2)),x,S(3),S(2)/S(7)*d*(c*d-b*e)*(d+e*x)**(S(7)/S(2))/e**S(3)-S(2)/S(9)*(S(2)*c*d-b*e)*(d+e*x)**(S(9)/S(2))/e**S(3)+S(2)/S(11)*c*(d+e*x)**(S(11)/S(2))/e**S(3)],
[S(1)/((d+e*x)**(S(3)/S(2))*(b*x+c*x**S(2))**S(3)),x,S(11),-S(3)/S(4)*(S(16)*c**S(2)*d**S(2)+S(12)*b*c*d*e+S(5)*b**S(2)*e**S(2))*arctanh(sqrt(d+e*x)/sqrt(d))/(b**S(5)*d**(S(7)/S(2)))+S(3)/S(4)*c**(S(7)/S(2))*(S(16)*c**S(2)*d**S(2)-S(44)*b*c*d*e+S(33)*b**S(2)*e**S(2))*arctanh(sqrt(c)*sqrt(d+e*x)/sqrt(c*d-b*e))/(b**S(5)*(c*d-b*e)**(S(7)/S(2)))+S(3)/S(4)*e*(c**S(2)*d**S(2)-b*c*d*e-b**S(2)*e**S(2))*(S(8)*c**S(2)*d**S(2)-S(8)*b*c*d*e+S(5)*b**S(2)*e**S(2))/(b**S(4)*d**S(3)*(c*d-b*e)**S(3)*sqrt(d+e*x))+S(1)/S(4)*c*(S(12)*c**S(2)*d**S(2)-S(5)*b*c*d*e-S(5)*b**S(2)*e**S(2))/(b**S(3)*d**S(2)*(c*d-b*e)*(b+c*x)**S(2)*sqrt(d+e*x))+(-S(1)/S(2))/(b*d*x**S(2)*(b+c*x)**S(2)*sqrt(d+e*x))+S(1)/S(4)*(S(8)*c*d+S(5)*b*e)/(b**S(2)*d**S(2)*x*(b+c*x)**S(2)*sqrt(d+e*x))+S(1)/S(4)*c*(S(2)*c*d-b*e)*(S(12)*c**S(2)*d**S(2)-S(12)*b*c*d*e-S(5)*b**S(2)*e**S(2))/(b**S(4)*d**S(2)*(c*d-b*e)**S(2)*(b+c*x)*sqrt(d+e*x))],
[S(1)/((d+e*x)**(S(5)/S(2))*(b*x+c*x**S(2))**S(3)),x,S(12),S(1)/S(12)*e*(S(72)*c**S(4)*d**S(4)-S(144)*b*c**S(3)*d**S(3)*e+S(27)*b**S(2)*c**S(2)*d**S(2)*e**S(2)+S(45)*b**S(3)*c*d*e**S(3)-S(35)*b**S(4)*e**S(4))/(b**S(4)*d**S(3)*(c*d-b*e)**S(3)*(d+e*x)**(S(3)/S(2)))+S(1)/S(4)*c*(S(12)*c**S(2)*d**S(2)-S(3)*b*c*d*e-S(7)*b**S(2)*e**S(2))/(b**S(3)*d**S(2)*(c*d-b*e)*(b+c*x)**S(2)*(d+e*x)**(S(3)/S(2)))+(-S(1)/S(2))/(b*d*x**S(2)*(b+c*x)**S(2)*(d+e*x)**(S(3)/S(2)))+S(1)/S(4)*(S(8)*c*d+S(7)*b*e)/(b**S(2)*d**S(2)*x*(b+c*x)**S(2)*(d+e*x)**(S(3)/S(2)))+S(1)/S(4)*c*(S(2)*c*d-b*e)*(S(12)*c**S(2)*d**S(2)-S(12)*b*c*d*e-S(7)*b**S(2)*e**S(2))/(b**S(4)*d**S(2)*(c*d-b*e)**S(2)*(b+c*x)*(d+e*x)**(S(3)/S(2)))-S(1)/S(4)*(S(48)*c**S(2)*d**S(2)+S(60)*b*c*d*e+S(35)*b**S(2)*e**S(2))*arctanh(sqrt(d+e*x)/sqrt(d))/(b**S(5)*d**(S(9)/S(2)))+S(1)/S(4)*c**(S(9)/S(2))*(S(48)*c**S(2)*d**S(2)-S(156)*b*c*d*e+S(143)*b**S(2)*e**S(2))*arctanh(sqrt(c)*sqrt(d+e*x)/sqrt(c*d-b*e))/(b**S(5)*(c*d-b*e)**(S(9)/S(2)))+S(1)/S(4)*e*(S(2)*c*d-b*e)*(S(12)*c**S(4)*d**S(4)-S(24)*b*c**S(3)*d**S(3)*e+S(2)*b**S(2)*c**S(2)*d**S(2)*e**S(2)+S(10)*b**S(3)*c*d*e**S(3)-S(35)*b**S(4)*e**S(4))/(b**S(4)*d**S(4)*(c*d-b*e)**S(4)*sqrt(d+e*x))],

# Integrands of the form (d+e x)**(m/S(2)) (b x+c x**S(2))**(p/S(2))

[(d+e*x)**(S(3)/S(2))*sqrt(b*x+c*x**S(2)),x,S(9),S(2)/S(7)*e*(b*x+c*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/c-S(2)/S(105)*(S(2)*c*d-b*e)*(S(3)*c**S(2)*d**S(2)-S(3)*b*c*d*e+S(8)*b**S(2)*e**S(2))*EllipticE(sqrt(c)*sqrt(x)/sqrt(-b),sqrt(b*e/(c*d)))*sqrt(-b)*sqrt(x)*sqrt(S(1)+c*x/b)*sqrt(d+e*x)/(c**(S(5)/S(2))*e**S(2)*sqrt(S(1)+e*x/d)*sqrt(b*x+c*x**S(2)))+S(4)/S(105)*d*(c*d-b*e)*(S(3)*c**S(2)*d**S(2)-S(3)*b*c*d*e+S(2)*b**S(2)*e**S(2))*EllipticF(sqrt(c)*sqrt(x)/sqrt(-b),sqrt(b*e/(c*d)))*sqrt(-b)*sqrt(x)*sqrt(S(1)+c*x/b)*sqrt(S(1)+e*x/d)/(c**(S(5)/S(2))*e**S(2)*sqrt(d+e*x)*sqrt(b*x+c*x**S(2)))+S(2)/S(105)*(S(3)*c**S(2)*d**S(2)+S(9)*b*c*d*e-S(4)*b**S(2)*e**S(2)+S(12)*c*e*(S(2)*c*d-b*e)*x)*sqrt(d+e*x)*sqrt(b*x+c*x**S(2))/(c**S(2)*e)],
[sqrt(d+e*x)*sqrt(b*x+c*x**S(2)),x,S(9),-S(4)/S(15)*(c**S(2)*d**S(2)-b*c*d*e+b**S(2)*e**S(2))*EllipticE(sqrt(c)*sqrt(x)/sqrt(-b),sqrt(b*e/(c*d)))*sqrt(-b)*sqrt(x)*sqrt(S(1)+c*x/b)*sqrt(d+e*x)/(c**(S(3)/S(2))*e**S(2)*sqrt(S(1)+e*x/d)*sqrt(b*x+c*x**S(2)))+S(2)/S(15)*d*(c*d-b*e)*(S(2)*c*d-b*e)*EllipticF(sqrt(c)*sqrt(x)/sqrt(-b),sqrt(b*e/(c*d)))*sqrt(-b)*sqrt(x)*sqrt(S(1)+c*x/b)*sqrt(S(1)+e*x/d)/(c**(S(3)/S(2))*e**S(2)*sqrt(d+e*x)*sqrt(b*x+c*x**S(2)))+S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(b*x+c*x**S(2))/e-S(2)/S(15)*(S(2)*c*d-b*e)*sqrt(d+e*x)*sqrt(b*x+c*x**S(2))/(c*e)],
[sqrt(d+e*x)/sqrt(-S(2)*x-S(3)*x**S(2)),x,S(4),-S(2)*EllipticE(sqrt(S(3)/S(2))*sqrt(-x),sqrt(S(2)/S(3)*e/d))*sqrt(d+e*x)/(sqrt(S(3))*sqrt(S(1)+e*x/d))],
[S(1)/(sqrt(d+e*x)*sqrt(-S(2)*x-S(3)*x**S(2))),x,S(4),-S(2)*EllipticF(sqrt(S(3)/S(2))*sqrt(-x),sqrt(S(2)/S(3)*e/d))*sqrt(S(1)+e*x/d)/(sqrt(S(3))*sqrt(d+e*x))],

#  Note: Integrands are equal. 
[sqrt(S(1)-x)/(sqrt(-x)*sqrt(S(1)+x)),x,S(1),-S(2)*EllipticE(sqrt(-x),I)],
[sqrt(S(1)-x)/sqrt(-x-x**S(2)),x,S(2),-S(2)*EllipticE(sqrt(-x),I)],

# Integrands of the form (d+e x)**m (b x+c x**S(2))**p when c d-b e=S(0)
[(d+e*x)**m/(c*d*x+c*e*x**S(2))**S(2),x,S(3),-e*(d+e*x)**(-S(1)+m)*hypergeom([S(2),-S(1)+m],[m],S(1)+e*x/d)/(c**S(2)*d**S(2)*(S(1)-m))],

# Integrands of the form (d+e x)**m (b x+c x**S(2))**p with m symbolic
[(d+e*x)**m*(b*x+c*x**S(2))**S(3),x,S(3),d**S(3)*(c*d-b*e)**S(3)*(d+e*x)**(S(1)+m)/(e**S(7)*(S(1)+m))-S(3)*d**S(2)*(c*d-b*e)**S(2)*(S(2)*c*d-b*e)*(d+e*x)**(S(2)+m)/(e**S(7)*(S(2)+m))+S(3)*d*(c*d-b*e)*(S(5)*c**S(2)*d**S(2)-S(5)*b*c*d*e+b**S(2)*e**S(2))*(d+e*x)**(S(3)+m)/(e**S(7)*(S(3)+m))-(S(2)*c*d-b*e)*(S(10)*c**S(2)*d**S(2)-S(10)*b*c*d*e+b**S(2)*e**S(2))*(d+e*x)**(S(4)+m)/(e**S(7)*(S(4)+m))+S(3)*c*(S(5)*c**S(2)*d**S(2)-S(5)*b*c*d*e+b**S(2)*e**S(2))*(d+e*x)**(S(5)+m)/(e**S(7)*(S(5)+m))-S(3)*c**S(2)*(S(2)*c*d-b*e)*(d+e*x)**(S(6)+m)/(e**S(7)*(S(6)+m))+c**S(3)*(d+e*x)**(S(7)+m)/(e**S(7)*(S(7)+m))],
[(d+e*x)**m*(b*x+c*x**S(2))**S(2),x,S(3),d**S(2)*(c*d-b*e)**S(2)*(d+e*x)**(S(1)+m)/(e**S(5)*(S(1)+m))-S(2)*d*(c*d-b*e)*(S(2)*c*d-b*e)*(d+e*x)**(S(2)+m)/(e**S(5)*(S(2)+m))+(S(6)*c**S(2)*d**S(2)-S(6)*b*c*d*e+b**S(2)*e**S(2))*(d+e*x)**(S(3)+m)/(e**S(5)*(S(3)+m))-S(2)*c*(S(2)*c*d-b*e)*(d+e*x)**(S(4)+m)/(e**S(5)*(S(4)+m))+c**S(2)*(d+e*x)**(S(5)+m)/(e**S(5)*(S(5)+m))],

# Integrands of the form (d+e x)**m (b x+c x**S(2))**p with p symbolic
[(d+e*x)**m*(b*x+c*x**S(2))**p,x,S(2),(d+e*x)**(S(1)+m)*(b*x+c*x**S(2))**p*AppellFS(1)(S(1)+m,-p,-p,S(2)+m,(d+e*x)/d,c*(d+e*x)/(c*d-b*e))/(e*(S(1)+m)*(-e*x/d)**p*(S(1)-c*(d+e*x)/(c*d-b*e))**p)],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when b=S(0)

# Integrands of the form (d+e x)**m (a+c x**S(2))**p

# p>S(0)
[(d+e*x)**S(4)*(a+c*x**S(2)),x,S(2),S(1)/S(5)*(c*d**S(2)+a*e**S(2))*(d+e*x)**S(5)/e**S(3)-S(1)/S(3)*c*d*(d+e*x)**S(6)/e**S(3)+S(1)/S(7)*c*(d+e*x)**S(7)/e**S(3)],
[(a+c*x**S(2))**S(3)/(d+e*x)**S(10),x,S(2),-S(1)/S(9)*(c*d**S(2)+a*e**S(2))**S(3)/(e**S(7)*(d+e*x)**S(9))+S(3)/S(4)*c*d*(c*d**S(2)+a*e**S(2))**S(2)/(e**S(7)*(d+e*x)**S(8))-S(3)/S(7)*c*(c*d**S(2)+a*e**S(2))*(S(5)*c*d**S(2)+a*e**S(2))/(e**S(7)*(d+e*x)**S(7))+S(2)/S(3)*c**S(2)*d*(S(5)*c*d**S(2)+S(3)*a*e**S(2))/(e**S(7)*(d+e*x)**S(6))-S(3)/S(5)*c**S(2)*(S(5)*c*d**S(2)+a*e**S(2))/(e**S(7)*(d+e*x)**S(5))+S(3)/S(2)*c**S(3)*d/(e**S(7)*(d+e*x)**S(4))-S(1)/S(3)*c**S(3)/(e**S(7)*(d+e*x)**S(3))],
[(c+d*x)*(a+b*x**S(2))**S(4),x,S(3),a**S(4)*c*x+S(4)/S(3)*a**S(3)*b*c*x**S(3)+S(6)/S(5)*a**S(2)*b**S(2)*c*x**S(5)+S(4)/S(7)*a*b**S(3)*c*x**S(7)+S(1)/S(9)*b**S(4)*c*x**S(9)+S(1)/S(10)*d*(a+b*x**S(2))**S(5)/b],

# p<S(0)
[(d+e*x)**S(4)/(a+c*x**S(2)),x,S(5),e**S(2)*(S(6)*c*d**S(2)-a*e**S(2))*x/c**S(2)+S(2)*d*e**S(3)*x**S(2)/c+S(1)/S(3)*e**S(4)*x**S(3)/c+S(2)*d*e*(c*d**S(2)-a*e**S(2))*log(a+c*x**S(2))/c**S(2)+(c**S(2)*d**S(4)-S(6)*a*c*d**S(2)*e**S(2)+a**S(2)*e**S(4))*arctan(x*sqrt(c)/sqrt(a))/(c**(S(5)/S(2))*sqrt(a))],
[(d+e*x)**S(3)/(a+c*x**S(2)),x,S(5),S(3)*d*e**S(2)*x/c+S(1)/S(2)*e**S(3)*x**S(2)/c+S(1)/S(2)*e*(S(3)*c*d**S(2)-a*e**S(2))*log(a+c*x**S(2))/c**S(2)+d*(c*d**S(2)-S(3)*a*e**S(2))*arctan(x*sqrt(c)/sqrt(a))/(c**(S(3)/S(2))*sqrt(a))],
[S(1)/((d+e*x)**S(2)*(a+c*x**S(2))**S(4)),x,S(8),S(1)/S(16)*e*(S(5)*c**S(3)*d**S(6)+S(23)*a*c**S(2)*d**S(4)*e**S(2)+S(47)*a**S(2)*c*d**S(2)*e**S(4)-S(35)*a**S(3)*e**S(6))/(a**S(3)*(c*d**S(2)+a*e**S(2))**S(4)*(d+e*x))+S(1)/S(6)*(a*e+c*d*x)/(a*(c*d**S(2)+a*e**S(2))*(d+e*x)*(a+c*x**S(2))**S(3))+S(1)/S(24)*(-a*e*(c*d**S(2)-S(7)*a*e**S(2))+c*d*(S(5)*c*d**S(2)+S(13)*a*e**S(2))*x)/(a**S(2)*(c*d**S(2)+a*e**S(2))**S(2)*(d+e*x)*(a+c*x**S(2))**S(2))+S(1)/S(48)*(-a*e*(S(5)*c*d**S(2)-S(7)*a*e**S(2))*(c*d**S(2)+S(5)*a*e**S(2))+S(3)*c*d*(S(5)*c**S(2)*d**S(4)+S(18)*a*c*d**S(2)*e**S(2)+S(29)*a**S(2)*e**S(4))*x)/(a**S(3)*(c*d**S(2)+a*e**S(2))**S(3)*(d+e*x)*(a+c*x**S(2)))+S(8)*c*d*e**S(7)*log(d+e*x)/(c*d**S(2)+a*e**S(2))**S(5)-S(4)*c*d*e**S(7)*log(a+c*x**S(2))/(c*d**S(2)+a*e**S(2))**S(5)+S(1)/S(16)*(S(5)*c**S(4)*d**S(8)+S(28)*a*c**S(3)*d**S(6)*e**S(2)+S(70)*a**S(2)*c**S(2)*d**S(4)*e**S(4)+S(140)*a**S(3)*c*d**S(2)*e**S(6)-S(35)*a**S(4)*e**S(8))*arctan(x*sqrt(c)/sqrt(a))*sqrt(c)/(a**(S(7)/S(2))*(c*d**S(2)+a*e**S(2))**S(5))],

# Integrands of the form (d+e x)**m (a+c x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**S(4)*sqrt(a+c*x**S(2)),x,S(6),S(3)/S(10)*d*e*(d+e*x)**S(2)*(a+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(6)*e*(d+e*x)**S(3)*(a+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(120)*e*(S(8)*d*(S(13)*c*d**S(2)-S(8)*a*e**S(2))+S(3)*e*(S(16)*c*d**S(2)-S(5)*a*e**S(2))*x)*(a+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(16)*a*(S(8)*c**S(2)*d**S(4)-S(12)*a*c*d**S(2)*e**S(2)+a**S(2)*e**S(4))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(5)/S(2))+S(1)/S(16)*(S(8)*c**S(2)*d**S(4)-S(12)*a*c*d**S(2)*e**S(2)+a**S(2)*e**S(4))*x*sqrt(a+c*x**S(2))/c**S(2)],
[(d+e*x)**S(3)*sqrt(a+c*x**S(2)),x,S(5),S(1)/S(5)*e*(d+e*x)**S(2)*(a+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(60)*e*(S(8)*(S(6)*c*d**S(2)-a*e**S(2))+S(21)*c*d*e*x)*(a+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(8)*a*d*(S(4)*c*d**S(2)-S(3)*a*e**S(2))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(3)/S(2))+S(1)/S(8)*d*(S(4)*c*d**S(2)-S(3)*a*e**S(2))*x*sqrt(a+c*x**S(2))/c],
[(d+e*x)**S(2)*sqrt(a+c*x**S(2)),x,S(5),S(5)/S(12)*d*e*(a+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(4)*e*(d+e*x)*(a+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(8)*a*(S(4)*c*d**S(2)-a*e**S(2))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(3)/S(2))+S(1)/S(8)*(S(4)*c*d**S(2)-a*e**S(2))*x*sqrt(a+c*x**S(2))/c],

# p<S(0)
[(d+e*x)**S(4)/sqrt(a+c*x**S(2)),x,S(5),S(1)/S(8)*(S(8)*c**S(2)*d**S(4)-S(24)*a*c*d**S(2)*e**S(2)+S(3)*a**S(2)*e**S(4))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(5)/S(2))+S(7)/S(12)*d*e*(d+e*x)**S(2)*sqrt(a+c*x**S(2))/c+S(1)/S(4)*e*(d+e*x)**S(3)*sqrt(a+c*x**S(2))/c+S(1)/S(24)*e*(S(4)*d*(S(19)*c*d**S(2)-S(16)*a*e**S(2))+e*(S(26)*c*d**S(2)-S(9)*a*e**S(2))*x)*sqrt(a+c*x**S(2))/c**S(2)],
[(d+e*x)**S(3)/sqrt(a+c*x**S(2)),x,S(4),S(1)/S(2)*d*(S(2)*c*d**S(2)-S(3)*a*e**S(2))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(3)/S(2))+S(1)/S(3)*e*(d+e*x)**S(2)*sqrt(a+c*x**S(2))/c+S(1)/S(6)*e*(S(4)*(S(4)*c*d**S(2)-a*e**S(2))+S(5)*c*d*e*x)*sqrt(a+c*x**S(2))/c**S(2)],
[(d+e*x)**S(2)/sqrt(a+c*x**S(2)),x,S(4),S(1)/S(2)*(S(2)*c*d**S(2)-a*e**S(2))*arctanh(x*sqrt(c)/sqrt(a+c*x**S(2)))/c**(S(3)/S(2))+S(3)/S(2)*d*e*sqrt(a+c*x**S(2))/c+S(1)/S(2)*e*(d+e*x)*sqrt(a+c*x**S(2))/c],

# Integrands of the form (d+e x)**(m/S(2)) (a+c x**S(2))**p

# p>S(0)
[(d+e*x)**(S(5)/S(2))*(a+c*x**S(2)),x,S(2),S(2)/S(7)*(c*d**S(2)+a*e**S(2))*(d+e*x)**(S(7)/S(2))/e**S(3)-S(4)/S(9)*c*d*(d+e*x)**(S(9)/S(2))/e**S(3)+S(2)/S(11)*c*(d+e*x)**(S(11)/S(2))/e**S(3)],
[(d+e*x)**(S(3)/S(2))*(a+c*x**S(2)),x,S(2),S(2)/S(5)*(c*d**S(2)+a*e**S(2))*(d+e*x)**(S(5)/S(2))/e**S(3)-S(4)/S(7)*c*d*(d+e*x)**(S(7)/S(2))/e**S(3)+S(2)/S(9)*c*(d+e*x)**(S(9)/S(2))/e**S(3)],
[(a+c*x**S(2))*sqrt(d+e*x),x,S(2),S(2)/S(3)*(c*d**S(2)+a*e**S(2))*(d+e*x)**(S(3)/S(2))/e**S(3)-S(4)/S(5)*c*d*(d+e*x)**(S(5)/S(2))/e**S(3)+S(2)/S(7)*c*(d+e*x)**(S(7)/S(2))/e**S(3)],

# p<S(0)
[(d+e*x)**(S(5)/S(2))/(a-c*x**S(2)),x,S(6),-S(2)/S(3)*e*(d+e*x)**(S(3)/S(2))/c-arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(-e*sqrt(a)+d*sqrt(c)))*(-e*sqrt(a)+d*sqrt(c))**(S(5)/S(2))/(c**(S(7)/S(4))*sqrt(a))+arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(e*sqrt(a)+d*sqrt(c)))*(e*sqrt(a)+d*sqrt(c))**(S(5)/S(2))/(c**(S(7)/S(4))*sqrt(a))-S(4)*d*e*sqrt(d+e*x)/c],
[(d+e*x)**(S(3)/S(2))/(a-c*x**S(2)),x,S(5),-arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(-e*sqrt(a)+d*sqrt(c)))*(-e*sqrt(a)+d*sqrt(c))**(S(3)/S(2))/(c**(S(5)/S(4))*sqrt(a))+arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(e*sqrt(a)+d*sqrt(c)))*(e*sqrt(a)+d*sqrt(c))**(S(3)/S(2))/(c**(S(5)/S(4))*sqrt(a))-S(2)*e*sqrt(d+e*x)/c],
[sqrt(d+e*x)/(a-c*x**S(2)),x,S(4),-arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(-e*sqrt(a)+d*sqrt(c)))*sqrt(-e*sqrt(a)+d*sqrt(c))/(c**(S(3)/S(4))*sqrt(a))+arctanh(c**(S(1)/S(4))*sqrt(d+e*x)/sqrt(e*sqrt(a)+d*sqrt(c)))*sqrt(e*sqrt(a)+d*sqrt(c))/(c**(S(3)/S(4))*sqrt(a))],

# Integrands of the form (d+e x)**(m/S(2)) (a+c x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**(S(3)/S(2))*sqrt(a+c*x**S(2)),x,S(7),S(2)/S(7)*e*(a+c*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/c+S(2)/S(105)*(S(3)*c*d**S(2)-S(5)*a*e**S(2)+S(24)*c*d*e*x)*sqrt(d+e*x)*sqrt(a+c*x**S(2))/(c*e)+S(4)/S(105)*d*(S(3)*c*d**S(2)-S(29)*a*e**S(2))*EllipticE(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(d+e*x)*sqrt(S(1)+c*x**S(2)/a)/(e**S(2)*sqrt(c)*sqrt(a+c*x**S(2))*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c))))-S(4)/S(105)*(S(3)*c*d**S(2)-S(5)*a*e**S(2))*(c*d**S(2)+a*e**S(2))*EllipticF(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(S(1)+c*x**S(2)/a)*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c)))/(c**(S(3)/S(2))*e**S(2)*sqrt(d+e*x)*sqrt(a+c*x**S(2)))],
[sqrt(d+e*x)*sqrt(a+c*x**S(2)),x,S(7),S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(a+c*x**S(2))/e-S(4)/S(15)*d*sqrt(d+e*x)*sqrt(a+c*x**S(2))/e+S(4)/S(15)*(c*d**S(2)-S(3)*a*e**S(2))*EllipticE(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(d+e*x)*sqrt(S(1)+c*x**S(2)/a)/(e**S(2)*sqrt(c)*sqrt(a+c*x**S(2))*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c))))-S(4)/S(15)*d*(c*d**S(2)+a*e**S(2))*EllipticF(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(S(1)+c*x**S(2)/a)*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c)))/(e**S(2)*sqrt(c)*sqrt(d+e*x)*sqrt(a+c*x**S(2)))],
[(a+c*x**S(2))**(S(5)/S(2))/(d+e*x)**(S(11)/S(2)),x,S(8),-S(4)/S(63)*c*(S(2)*d*(S(4)*c*d**S(2)+a*e**S(2))+e*(S(13)*c*d**S(2)+S(7)*a*e**S(2))*x)*(a+c*x**S(2))**(S(3)/S(2))/(e**S(3)*(c*d**S(2)+a*e**S(2))*(d+e*x)**(S(7)/S(2)))-S(2)/S(9)*(a+c*x**S(2))**(S(5)/S(2))/(e*(d+e*x)**(S(9)/S(2)))-S(8)/S(63)*c**S(2)*(d*(S(32)*c**S(2)*d**S(4)+S(49)*a*c*d**S(2)*e**S(2)+S(9)*a**S(2)*e**S(4))+e*(S(40)*c**S(2)*d**S(4)+S(69)*a*c*d**S(2)*e**S(2)+S(21)*a**S(2)*e**S(4))*x)*sqrt(a+c*x**S(2))/(e**S(5)*(c*d**S(2)+a*e**S(2))**S(2)*(d+e*x)**(S(3)/S(2)))-S(16)/S(63)*c**(S(5)/S(2))*(S(32)*c**S(2)*d**S(4)+S(57)*a*c*d**S(2)*e**S(2)+S(21)*a**S(2)*e**S(4))*EllipticE(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(d+e*x)*sqrt(S(1)+c*x**S(2)/a)/(e**S(6)*(c*d**S(2)+a*e**S(2))**S(2)*sqrt(a+c*x**S(2))*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c))))+S(16)/S(63)*c**(S(5)/S(2))*d*(S(32)*c*d**S(2)+S(33)*a*e**S(2))*EllipticF(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(S(1)+c*x**S(2)/a)*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c)))/(e**S(6)*(c*d**S(2)+a*e**S(2))*sqrt(d+e*x)*sqrt(a+c*x**S(2)))],

# p<S(0)
[(d+e*x)**(S(7)/S(2))/sqrt(a+c*x**S(2)),x,S(8),S(24)/S(35)*d*e*(d+e*x)**(S(3)/S(2))*sqrt(a+c*x**S(2))/c+S(2)/S(7)*e*(d+e*x)**(S(5)/S(2))*sqrt(a+c*x**S(2))/c+S(2)/S(105)*e*(S(71)*c*d**S(2)-S(25)*a*e**S(2))*sqrt(d+e*x)*sqrt(a+c*x**S(2))/c**S(2)-S(32)/S(105)*d*(S(11)*c*d**S(2)-S(13)*a*e**S(2))*EllipticE(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(d+e*x)*sqrt(S(1)+c*x**S(2)/a)/(c**(S(3)/S(2))*sqrt(a+c*x**S(2))*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c))))+S(2)/S(105)*(S(71)*c*d**S(2)-S(25)*a*e**S(2))*(c*d**S(2)+a*e**S(2))*EllipticF(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(S(1)+c*x**S(2)/a)*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c)))/(c**(S(5)/S(2))*sqrt(d+e*x)*sqrt(a+c*x**S(2)))],
[(d+e*x)**(S(5)/S(2))/sqrt(a+c*x**S(2)),x,S(7),S(2)/S(5)*e*(d+e*x)**(S(3)/S(2))*sqrt(a+c*x**S(2))/c+S(16)/S(15)*d*e*sqrt(d+e*x)*sqrt(a+c*x**S(2))/c-S(2)/S(15)*(S(23)*c*d**S(2)-S(9)*a*e**S(2))*EllipticE(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(d+e*x)*sqrt(S(1)+c*x**S(2)/a)/(c**(S(3)/S(2))*sqrt(a+c*x**S(2))*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c))))+S(16)/S(15)*d*(c*d**S(2)+a*e**S(2))*EllipticF(sqrt(S(1)-x*sqrt(c)/sqrt(-a))/sqrt(S(2)),sqrt(-S(2)*a*e/(-a*e+d*sqrt(-a)*sqrt(c))))*sqrt(-a)*sqrt(S(1)+c*x**S(2)/a)*sqrt((d+e*x)*sqrt(c)/(e*sqrt(-a)+d*sqrt(c)))/(c**(S(3)/S(2))*sqrt(d+e*x)*sqrt(a+c*x**S(2)))],

# Integrands of the form (d+e x)**m (a+c x**S(2))**(p/S(3))
[S(1)/((d+e*x)*(d**S(2)+S(3)*e**S(2)*x**S(2))**(S(1)/S(3))),x,S(1),-S(1)/S(2)*log(d+e*x)/(S(2)**(S(2)/S(3))*d**(S(2)/S(3))*e)+S(1)/S(2)*log(S(3)*d*e**S(2)-S(3)*e**S(3)*x-S(3)*S(2)**(S(1)/S(3))*d**(S(1)/S(3))*e**S(2)*(d**S(2)+S(3)*e**S(2)*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*d**(S(2)/S(3))*e)-arctan(S(1)/sqrt(S(3))+S(2)**(S(2)/S(3))*(d-e*x)/(d**(S(1)/S(3))*(d**S(2)+S(3)*e**S(2)*x**S(2))**(S(1)/S(3))*sqrt(S(3))))/(S(2)**(S(2)/S(3))*d**(S(2)/S(3))*e*sqrt(S(3)))],
[(S(2)+S(3)*x)**S(3)/(S(4)+S(27)*x**S(2))**(S(1)/S(3)),x,S(6),S(1)/S(30)*(S(2)+S(3)*x)**S(2)*(S(4)+S(27)*x**S(2))**(S(2)/S(3))+S(4)/S(35)*(S(7)+S(4)*x)*(S(4)+S(27)*x**S(2))**(S(2)/S(3))-S(96)/S(7)*x/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3))))-S(32)/S(63)*S(2)**(S(5)/S(6))*(S(2)**(S(2)/S(3))-(S(4)+S(27)*x**S(2))**(S(1)/S(3)))*EllipticF((-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)+sqrt(S(3))))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3)))),sqrt(-S(7)+S(4)*sqrt(S(3))))*sqrt((S(2)*S(2)**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(4)+S(27)*x**S(2))**(S(1)/S(3))+(S(4)+S(27)*x**S(2))**(S(2)/S(3)))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3))))**S(2))/(S(3)**(S(1)/S(4))*x*sqrt((-S(2)**(S(2)/S(3))+(S(4)+S(27)*x**S(2))**(S(1)/S(3)))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3))))**S(2)))+S(16)/S(21)*S(2)**(S(1)/S(3))*(S(2)**(S(2)/S(3))-(S(4)+S(27)*x**S(2))**(S(1)/S(3)))*EllipticE((-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)+sqrt(S(3))))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3)))),sqrt(-S(7)+S(4)*sqrt(S(3))))*sqrt((S(2)*S(2)**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(4)+S(27)*x**S(2))**(S(1)/S(3))+(S(4)+S(27)*x**S(2))**(S(2)/S(3)))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3))))**S(2))*sqrt(S(2)+sqrt(S(3)))/(S(3)**(S(3)/S(4))*x*sqrt((-S(2)**(S(2)/S(3))+(S(4)+S(27)*x**S(2))**(S(1)/S(3)))/(-(S(4)+S(27)*x**S(2))**(S(1)/S(3))+S(2)**(S(2)/S(3))*(S(1)-sqrt(S(3))))**S(2)))],

# Integrands of the form (d+e x)**m (a+c x**S(2))**(p/S(4))
[S(1)/((a+b*x)*(c+d*x**S(2))**(S(1)/S(4))),x,S(11),arctan((c+d*x**S(2))**(S(1)/S(4))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(1)/S(4)))/((b**S(2)*c+a**S(2)*d)**(S(1)/S(4))*sqrt(b))-arctanh((c+d*x**S(2))**(S(1)/S(4))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(1)/S(4)))/((b**S(2)*c+a**S(2)*d)**(S(1)/S(4))*sqrt(b))-a*c**(S(1)/S(4))*EllipticPi((c+d*x**S(2))**(S(1)/S(4))/c**(S(1)/S(4)),-b*sqrt(c)/sqrt(b**S(2)*c+a**S(2)*d),I)*sqrt(-d*x**S(2)/c)/(b*x*sqrt(b**S(2)*c+a**S(2)*d))+a*c**(S(1)/S(4))*EllipticPi((c+d*x**S(2))**(S(1)/S(4))/c**(S(1)/S(4)),b*sqrt(c)/sqrt(b**S(2)*c+a**S(2)*d),I)*sqrt(-d*x**S(2)/c)/(b*x*sqrt(b**S(2)*c+a**S(2)*d))],
[S(1)/((a+b*x)*(c+d*x**S(2))**(S(3)/S(4))),x,S(11),-arctan((c+d*x**S(2))**(S(1)/S(4))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(1)/S(4)))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(3)/S(4))-arctanh((c+d*x**S(2))**(S(1)/S(4))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(1)/S(4)))*sqrt(b)/(b**S(2)*c+a**S(2)*d)**(S(3)/S(4))+a*c**(S(1)/S(4))*EllipticPi((c+d*x**S(2))**(S(1)/S(4))/c**(S(1)/S(4)),-b*sqrt(c)/sqrt(b**S(2)*c+a**S(2)*d),I)*sqrt(-d*x**S(2)/c)/((b**S(2)*c+a**S(2)*d)*x)+a*c**(S(1)/S(4))*EllipticPi((c+d*x**S(2))**(S(1)/S(4))/c**(S(1)/S(4)),b*sqrt(c)/sqrt(b**S(2)*c+a**S(2)*d),I)*sqrt(-d*x**S(2)/c)/((b**S(2)*c+a**S(2)*d)*x)],
[S(1)/((d+e*x)**(S(3)/S(2))*(a+c*x**S(2))**(S(1)/S(4))),x,S(1),-S(2)*hypergeom([-S(1)/S(2),S(1)/S(4)],[S(1)/S(2)],S(2)*(d+e*x)*sqrt(-a)*sqrt(c)/((-e*sqrt(-a)+d*sqrt(c))*(sqrt(-a)-x*sqrt(c))))*(sqrt(-a)-x*sqrt(c))*(-(e*sqrt(-a)+d*sqrt(c))*(sqrt(-a)+x*sqrt(c))/((-e*sqrt(-a)+d*sqrt(c))*(sqrt(-a)-x*sqrt(c))))**(S(1)/S(4))/((a+c*x**S(2))**(S(1)/S(4))*(e*sqrt(-a)+d*sqrt(c))*sqrt(d+e*x))],

# Integrands of the form (d+e x)**m (a+c x**S(2))**p with m symbolic
[(d+e*x)**m*(a+c*x**S(2))**S(3),x,S(2),(c*d**S(2)+a*e**S(2))**S(3)*(d+e*x)**(S(1)+m)/(e**S(7)*(S(1)+m))-S(6)*c*d*(c*d**S(2)+a*e**S(2))**S(2)*(d+e*x)**(S(2)+m)/(e**S(7)*(S(2)+m))+S(3)*c*(c*d**S(2)+a*e**S(2))*(S(5)*c*d**S(2)+a*e**S(2))*(d+e*x)**(S(3)+m)/(e**S(7)*(S(3)+m))-S(4)*c**S(2)*d*(S(5)*c*d**S(2)+S(3)*a*e**S(2))*(d+e*x)**(S(4)+m)/(e**S(7)*(S(4)+m))+S(3)*c**S(2)*(S(5)*c*d**S(2)+a*e**S(2))*(d+e*x)**(S(5)+m)/(e**S(7)*(S(5)+m))-S(6)*c**S(3)*d*(d+e*x)**(S(6)+m)/(e**S(7)*(S(6)+m))+c**S(3)*(d+e*x)**(S(7)+m)/(e**S(7)*(S(7)+m))],
[(d+e*x)**m*(a+c*x**S(2))**S(2),x,S(2),(c*d**S(2)+a*e**S(2))**S(2)*(d+e*x)**(S(1)+m)/(e**S(5)*(S(1)+m))-S(4)*c*d*(c*d**S(2)+a*e**S(2))*(d+e*x)**(S(2)+m)/(e**S(5)*(S(2)+m))+S(2)*c*(S(3)*c*d**S(2)+a*e**S(2))*(d+e*x)**(S(3)+m)/(e**S(5)*(S(3)+m))-S(4)*c**S(2)*d*(d+e*x)**(S(4)+m)/(e**S(5)*(S(4)+m))+c**S(2)*(d+e*x)**(S(5)+m)/(e**S(5)*(S(5)+m))],

# Integrands of the form (d+e x)**m (a+c x**S(2))**p with p symbolic
[(d+e*x)**m*(a+c*x**S(2))**p,x,S(2),(d+e*x)**(S(1)+m)*(a+c*x**S(2))**p*AppellFS(1)(S(1)+m,-p,-p,S(2)+m,(d+e*x)/(d-e*sqrt(-a)/sqrt(c)),(d+e*x)/(d+e*sqrt(-a)/sqrt(c)))/(e*(S(1)+m)*(S(1)+(-d-e*x)/(d-e*sqrt(-a)/sqrt(c)))**p*(S(1)+(-d-e*x)/(d+e*sqrt(-a)/sqrt(c)))**p)],
[(d+e*x)**S(3)*(a+c*x**S(2))**p,x,S(4),S(1)/S(2)*e*(d+e*x)**S(2)*(a+c*x**S(2))**(S(1)+p)/(c*(S(2)+p))-S(1)/S(2)*e*((S(3)+S(2)*p)*(a*e**S(2)-c*d**S(2)*(S(5)+S(2)*p))-S(2)*c*d*e*(S(1)+p)*(S(3)+p)*x)*(a+c*x**S(2))**(S(1)+p)/(c**S(2)*(S(2)+p)*(S(3)+S(5)*p+S(2)*p**S(2)))-d*(S(3)*a*e**S(2)-c*d**S(2)*(S(3)+S(2)*p))*x*(a+c*x**S(2))**p*hypergeom([S(1)/S(2),-p],[S(3)/S(2)],-c*x**S(2)/a)/(c*(S(3)+S(2)*p)*(S(1)+c*x**S(2)/a)**p)],
[(d+e*x)**S(2)*(a+c*x**S(2))**p,x,S(4),d*e*(a+c*x**S(2))**(S(1)+p)/(c*(S(1)+p))+e**S(2)*x*(a+c*x**S(2))**(S(1)+p)/(c*(S(3)+S(2)*p))-(a*e**S(2)-c*d**S(2)*(S(3)+S(2)*p))*x*(a+c*x**S(2))**p*hypergeom([S(1)/S(2),-p],[S(3)/S(2)],-c*x**S(2)/a)/(c*(S(3)+S(2)*p)*(S(1)+c*x**S(2)/a)**p),d*e*(S(2)+p)*(a+c*x**S(2))**(S(1)+p)/(c*(S(3)+S(5)*p+S(2)*p**S(2)))+e*(d+e*x)*(a+c*x**S(2))**(S(1)+p)/(c*(S(3)+S(2)*p))-(a*e**S(2)-c*d**S(2)*(S(3)+S(2)*p))*x*(a+c*x**S(2))**p*hypergeom([S(1)/S(2),-p],[S(3)/S(2)],-c*x**S(2)/a)/(c*(S(3)+S(2)*p)*(S(1)+c*x**S(2)/a)**p)],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when b=S(0) and c d**S(2)+a e**S(2)=S(0)

# Integrands of the form (d+e x)**m (d**S(2) - e**S(2) x**S(2))**p

# p>S(0)

# p<S(0)
[(a+b*x)**S(6)/(a**S(2)-b**S(2)*x**S(2)),x,S(3),-S(16)*a**S(4)*x-S(4)*a**S(3)*(a+b*x)**S(2)/b-S(4)/S(3)*a**S(2)*(a+b*x)**S(3)/b-S(1)/S(2)*a*(a+b*x)**S(4)/b-S(1)/S(5)*(a+b*x)**S(5)/b-S(32)*a**S(5)*log(a-b*x)/b],
[(a+b*x)**S(5)/(a**S(2)-b**S(2)*x**S(2)),x,S(3),-S(8)*a**S(3)*x-S(2)*a**S(2)*(a+b*x)**S(2)/b-S(2)/S(3)*a*(a+b*x)**S(3)/b-S(1)/S(4)*(a+b*x)**S(4)/b-S(16)*a**S(4)*log(a-b*x)/b],
[(a+b*x)**S(4)/(a**S(2)-b**S(2)*x**S(2)),x,S(3),-S(4)*a**S(2)*x-a*(a+b*x)**S(2)/b-S(1)/S(3)*(a+b*x)**S(3)/b-S(8)*a**S(3)*log(a-b*x)/b],

# Integrands of the form (d+e x)**m (d**S(2) - e**S(2) x**S(2))**(p/S(2))

# p>S(0)
[(a+b*x)**S(4)*sqrt(a**S(2)-b**S(2)*x**S(2)),x,S(7),-S(7)/S(8)*a**S(3)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(21)/S(40)*a**S(2)*(a+b*x)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(3)/S(10)*a*(a+b*x)**S(2)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(1)/S(6)*(a+b*x)**S(3)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b+S(21)/S(16)*a**S(6)*arctan(b*x/sqrt(a**S(2)-b**S(2)*x**S(2)))/b+S(21)/S(16)*a**S(4)*x*sqrt(a**S(2)-b**S(2)*x**S(2))],
[(a+b*x)**S(3)*sqrt(a**S(2)-b**S(2)*x**S(2)),x,S(6),-S(7)/S(12)*a**S(2)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(7)/S(20)*a*(a+b*x)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(1)/S(5)*(a+b*x)**S(2)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b+S(7)/S(8)*a**S(5)*arctan(b*x/sqrt(a**S(2)-b**S(2)*x**S(2)))/b+S(7)/S(8)*a**S(3)*x*sqrt(a**S(2)-b**S(2)*x**S(2))],
[(a+b*x)**S(2)*sqrt(a**S(2)-b**S(2)*x**S(2)),x,S(5),-S(5)/S(12)*a*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b-S(1)/S(4)*(a+b*x)*(a**S(2)-b**S(2)*x**S(2))**(S(3)/S(2))/b+S(5)/S(8)*a**S(4)*arctan(b*x/sqrt(a**S(2)-b**S(2)*x**S(2)))/b+S(5)/S(8)*a**S(2)*x*sqrt(a**S(2)-b**S(2)*x**S(2))],

# p<S(0)
[(d+e*x)**S(5)/sqrt(d**S(2)-e**S(2)*x**S(2)),x,S(7),S(63)/S(8)*d**S(5)*arctan(e*x/sqrt(d**S(2)-e**S(2)*x**S(2)))/e-S(63)/S(8)*d**S(4)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(21)/S(8)*d**S(3)*(d+e*x)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(21)/S(20)*d**S(2)*(d+e*x)**S(2)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(9)/S(20)*d*(d+e*x)**S(3)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(1)/S(5)*(d+e*x)**S(4)*sqrt(d**S(2)-e**S(2)*x**S(2))/e],
[(d+e*x)**S(4)/sqrt(d**S(2)-e**S(2)*x**S(2)),x,S(6),S(35)/S(8)*d**S(4)*arctan(e*x/sqrt(d**S(2)-e**S(2)*x**S(2)))/e-S(35)/S(8)*d**S(3)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(35)/S(24)*d**S(2)*(d+e*x)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(7)/S(12)*d*(d+e*x)**S(2)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(1)/S(4)*(d+e*x)**S(3)*sqrt(d**S(2)-e**S(2)*x**S(2))/e],
[(d+e*x)**S(3)/sqrt(d**S(2)-e**S(2)*x**S(2)),x,S(5),S(5)/S(2)*d**S(3)*arctan(e*x/sqrt(d**S(2)-e**S(2)*x**S(2)))/e-S(5)/S(2)*d**S(2)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(5)/S(6)*d*(d+e*x)*sqrt(d**S(2)-e**S(2)*x**S(2))/e-S(1)/S(3)*(d+e*x)**S(2)*sqrt(d**S(2)-e**S(2)*x**S(2))/e],

# Integrands of the form (d+e x)**(m/S(2)) (c d**S(2) - c e**S(2) x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**(S(5)/S(2))*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(4),-S(256)/S(315)*d**S(3)*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*(d+e*x)**(S(3)/S(2)))-S(2)/S(9)*(d+e*x)**(S(3)/S(2))*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e)-S(64)/S(105)*d**S(2)*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*sqrt(d+e*x))-S(8)/S(21)*d*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/(c*e)],
[(d+e*x)**(S(3)/S(2))*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),-S(64)/S(105)*d**S(2)*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*(d+e*x)**(S(3)/S(2)))-S(16)/S(35)*d*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*sqrt(d+e*x))-S(2)/S(7)*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/(c*e)],
[(d+e*x)**(S(1)/S(2))*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),-S(8)/S(15)*d*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*(d+e*x)**(S(3)/S(2)))-S(2)/S(5)*(c*d**S(2)-c*e**S(2)*x**S(2))**(S(3)/S(2))/(c*e*sqrt(d+e*x))],

# p<S(0)
[(d+e*x)**(S(7)/S(2))/(c*d**S(2)-c*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(4),-S(24)/S(35)*d*(d+e*x)**(S(3)/S(2))*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e)-S(2)/S(7)*(d+e*x)**(S(5)/S(2))*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e)-S(256)/S(35)*d**S(3)*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e*sqrt(d+e*x))-S(64)/S(35)*d**S(2)*sqrt(d+e*x)*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e)],
[(d+e*x)**(S(5)/S(2))/(c*d**S(2)-c*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),-S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e)-S(64)/S(15)*d**S(2)*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e*sqrt(d+e*x))-S(16)/S(15)*d*sqrt(d+e*x)*sqrt(c*d**S(2)-c*e**S(2)*x**S(2))/(c*e)],

# Integrands of the form (d+e x)**(m/S(2)) (c d**S(2) - c e**S(2) x**S(2))**(p/S(2)) when d>S(0) and c>S(0)

# p>S(0)
[(S(2)+e*x)**(S(5)/S(2))*(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),-S(128)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))+S(2)/S(3)*(S(2)-e*x)**(S(9)/S(2))/(e*sqrt(S(3)))+S(96)/S(5)*(S(2)-e*x)**(S(5)/S(2))*sqrt(S(3))/e-S(24)/S(7)*(S(2)-e*x)**(S(7)/S(2))*sqrt(S(3))/e],
[(S(2)+e*x)**(S(3)/S(2))*(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),-S(32)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))+S(16)/S(5)*(S(2)-e*x)**(S(5)/S(2))*sqrt(S(3))/e-S(2)/S(7)*(S(2)-e*x)**(S(7)/S(2))*sqrt(S(3))/e],
[(S(2)+e*x)**(S(1)/S(2))*(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),-S(8)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))+S(2)/S(5)*(S(2)-e*x)**(S(5)/S(2))*sqrt(S(3))/e],

# p<S(0)
[(S(2)+e*x)**(S(7)/S(2))/(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),S(32)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))+S(2)/S(7)*(S(2)-e*x)**(S(7)/S(2))/(e*sqrt(S(3)))-S(8)/S(5)*(S(2)-e*x)**(S(5)/S(2))*sqrt(S(3))/e-S(128)*sqrt(S(2)-e*x)/(e*sqrt(S(3)))],
[(S(2)+e*x)**(S(5)/S(2))/(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),S(16)/S(3)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))-S(2)/S(5)*(S(2)-e*x)**(S(5)/S(2))/(e*sqrt(S(3)))-S(32)*sqrt(S(2)-e*x)/(e*sqrt(S(3)))],
[(S(2)+e*x)**(S(3)/S(2))/(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(2)),x,S(3),S(2)/S(3)*(S(2)-e*x)**(S(3)/S(2))/(e*sqrt(S(3)))-S(8)*sqrt(S(2)-e*x)/(e*sqrt(S(3)))],

#  The following pairs of integrands are equal: 
[S(1)/((S(1)+a*x)*sqrt(S(1)-a*x)),x,S(2),-arctanh(sqrt(S(1)-a*x)/sqrt(S(2)))*sqrt(S(2))/a],
[S(1)/(sqrt(S(1)+a*x)*sqrt(S(1)-a**S(2)*x**S(2))),x,S(3),-arctanh(sqrt(S(1)-a*x)/sqrt(S(2)))*sqrt(S(2))/a],

# Integrands of the form (d+e x)**(m/S(2)) (c d**S(2) - c e**S(2) x**S(2))**(p/S(4)) when d>S(0) and c>S(0)

# p>S(0)
[(S(2)+e*x)**(S(1)/S(2))*(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(4)),x,S(15),S(3)/S(2)*S(3)**(S(1)/S(4))*(S(2)-e*x)**(S(1)/S(4))*(S(2)+e*x)**(S(3)/S(4))/e-S(1)/S(2)*S(3)**(S(1)/S(4))*(S(2)-e*x)**(S(5)/S(4))*(S(2)+e*x)**(S(3)/S(4))/e+S(3)*S(3)**(S(1)/S(4))*arctan(S(1)-(S(2)-e*x)**(S(1)/S(4))*sqrt(S(2))/(S(2)+e*x)**(S(1)/S(4)))/(e*sqrt(S(2)))-S(3)*S(3)**(S(1)/S(4))*arctan(S(1)+(S(2)-e*x)**(S(1)/S(4))*sqrt(S(2))/(S(2)+e*x)**(S(1)/S(4)))/(e*sqrt(S(2)))+S(3)/S(2)*S(3)**(S(1)/S(4))*log(sqrt(S(3))-(S(2)-e*x)**(S(1)/S(4))*sqrt(S(6))/(S(2)+e*x)**(S(1)/S(4))+sqrt(S(3))*sqrt(S(2)-e*x)/sqrt(S(2)+e*x))/(e*sqrt(S(2)))-S(3)/S(2)*S(3)**(S(1)/S(4))*log(sqrt(S(3))+(S(2)-e*x)**(S(1)/S(4))*sqrt(S(6))/(S(2)+e*x)**(S(1)/S(4))+sqrt(S(3))*sqrt(S(2)-e*x)/sqrt(S(2)+e*x))/(e*sqrt(S(2)))],
[(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(4))/(S(2)+e*x)**(S(9)/S(2)),x,S(3),-S(1)/S(13)*S(3)**(S(1)/S(4))*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(e*(S(2)+e*x)**(S(9)/S(2)))-S(2)/S(39)*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(S(3)**(S(3)/S(4))*e*(S(2)+e*x)**(S(7)/S(2)))-S(2)/S(195)*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(S(3)**(S(3)/S(4))*e*(S(2)+e*x)**(S(5)/S(2)))],
[(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(4))/(S(2)+e*x)**(S(11)/S(2)),x,S(4),-S(1)/S(17)*S(3)**(S(1)/S(4))*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(e*(S(2)+e*x)**(S(11)/S(2)))-S(3)/S(221)*S(3)**(S(1)/S(4))*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(e*(S(2)+e*x)**(S(9)/S(2)))-S(2)/S(221)*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(S(3)**(S(3)/S(4))*e*(S(2)+e*x)**(S(7)/S(2)))-S(2)/S(1105)*(S(4)-e**S(2)*x**S(2))**(S(5)/S(4))/(S(3)**(S(3)/S(4))*e*(S(2)+e*x)**(S(5)/S(2)))],

# p<S(0)
[(S(2)+e*x)**(S(5)/S(2))/(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(4)),x,S(16),-S(5)/S(2)*S(3)**(S(3)/S(4))*(S(2)-e*x)**(S(3)/S(4))*(S(2)+e*x)**(S(1)/S(4))/e-S(1)/S(2)*S(3)**(S(3)/S(4))*(S(2)-e*x)**(S(3)/S(4))*(S(2)+e*x)**(S(5)/S(4))/e-S(1)/S(3)*(S(2)-e*x)**(S(3)/S(4))*(S(2)+e*x)**(S(9)/S(4))/(S(3)**(S(1)/S(4))*e)+S(5)*S(3)**(S(3)/S(4))*arctan(S(1)-(S(2)-e*x)**(S(1)/S(4))*sqrt(S(2))/(S(2)+e*x)**(S(1)/S(4)))/(e*sqrt(S(2)))-S(5)*S(3)**(S(3)/S(4))*arctan(S(1)+(S(2)-e*x)**(S(1)/S(4))*sqrt(S(2))/(S(2)+e*x)**(S(1)/S(4)))/(e*sqrt(S(2)))-S(5)/S(2)*S(3)**(S(3)/S(4))*log(sqrt(S(3))-(S(2)-e*x)**(S(1)/S(4))*sqrt(S(6))/(S(2)+e*x)**(S(1)/S(4))+sqrt(S(3))*sqrt(S(2)-e*x)/sqrt(S(2)+e*x))/(e*sqrt(S(2)))+S(5)/S(2)*S(3)**(S(3)/S(4))*log(sqrt(S(3))+(S(2)-e*x)**(S(1)/S(4))*sqrt(S(6))/(S(2)+e*x)**(S(1)/S(4))+sqrt(S(3))*sqrt(S(2)-e*x)/sqrt(S(2)+e*x))/(e*sqrt(S(2)))],
[S(1)/((S(2)+e*x)**(S(9)/S(2))*(S(12)-S(3)*e**S(2)*x**S(2))**(S(1)/S(4))),x,S(4),-S(1)/S(15)*(S(4)-e**S(2)*x**S(2))**(S(3)/S(4))/(S(3)**(S(1)/S(4))*e*(S(2)+e*x)**(S(9)/S(2)))-S(1)/S(55)*(S(4)-e**S(2)*x**S(2))**(S(3)/S(4))/(S(3)**(S(1)/S(4))*e*(S(2)+e*x)**(S(7)/S(2)))-S(2)/S(385)*(S(4)-e**S(2)*x**S(2))**(S(3)/S(4))/(S(3)**(S(1)/S(4))*e*(S(2)+e*x)**(S(5)/S(2)))-S(2)/S(1155)*(S(4)-e**S(2)*x**S(2))**(S(3)/S(4))/(S(3)**(S(1)/S(4))*e*(S(2)+e*x)**(S(3)/S(2)))],

# Integrands of the form (d+e x)**m (d**S(2) - e**S(2) x**S(2))**p with m symbolic
[(a+b*x)**m*(a**S(2)-b**S(2)*x**S(2))**S(3),x,S(3),S(8)*a**S(3)*(a+b*x)**(S(4)+m)/(b*(S(4)+m))-S(12)*a**S(2)*(a+b*x)**(S(5)+m)/(b*(S(5)+m))+S(6)*a*(a+b*x)**(S(6)+m)/(b*(S(6)+m))-(a+b*x)**(S(7)+m)/(b*(S(7)+m))],
[(d+e*x)**m/(d**S(2)-e**S(2)*x**S(2))**(S(7)/S(2)),x,S(3),-(d-e*x)*(d+e*x)**(S(1)+m)*hypergeom([S(1),-S(5)+m],[-S(3)/S(2)+m],S(1)/S(2)*(d+e*x)/d)/(d*e*(S(5)-S(2)*m)*(d**S(2)-e**S(2)*x**S(2))**(S(7)/S(2))),S(1)/S(5)*S(2)**(-S(5)/S(2)+m)*(d+e*x)**m*((d+e*x)/d)**(S(1)/S(2)-m)*hypergeom([-S(5)/S(2),S(7)/S(2)-m],[-S(3)/S(2)],S(1)/S(2)*(d-e*x)/d)/(d**S(3)*e*(d-e*x)**S(2)*sqrt(d**S(2)-e**S(2)*x**S(2)))],

# Integrands of the form (d+e x)**m (d**S(2) - e**S(2) x**S(2))**p with p symbolic
[(a+b*x)**m*(a**S(2)-b**S(2)*x**S(2))**p,x,S(3),S(1)/S(2)*(a-b*x)*(a+b*x)**(S(1)+m)*(a**S(2)-b**S(2)*x**S(2))**p*hypergeom([S(1),S(2)+m+S(2)*p],[S(2)+m+p],S(1)/S(2)*(a+b*x)/a)/(a*b*(S(1)+m+p)),-S(2)**(m+p)*(a-b*x)*(a+b*x)**m*((a+b*x)/a)**(-m-p)*(a**S(2)-b**S(2)*x**S(2))**p*hypergeom([-m-p,S(1)+p],[S(2)+p],S(1)/S(2)*(a-b*x)/a)/(b*(S(1)+p))],
[(a+b*x)**S(3)*(a**S(2)-b**S(2)*x**S(2))**p,x,S(3),S(1)/S(2)*(a-b*x)*(a+b*x)**S(4)*(a**S(2)-b**S(2)*x**S(2))**p*hypergeom([S(1),S(5)+S(2)*p],[S(5)+p],S(1)/S(2)*(a+b*x)/a)/(a*b*(S(4)+p)),-S(2)**(S(3)+p)*a**S(3)*(a-b*x)*(a**S(2)-b**S(2)*x**S(2))**p*hypergeom([-S(3)-p,S(1)+p],[S(2)+p],S(1)/S(2)*(a-b*x)/a)/(b*(S(1)+p)*((a+b*x)/a)**p)],
[(a**S(2)-b**S(2)*x**S(2))**p/(a+b*x)**(S(3)/S(2)),x,S(3),-S(2)**(-S(3)/S(2)+p)*(a-b*x)*((a+b*x)/a)**(S(1)/S(2)-p)*(a**S(2)-b**S(2)*x**S(2))**p*hypergeom([S(3)/S(2)-p,S(1)+p],[S(2)+p],S(1)/S(2)*(a-b*x)/a)/(a*b*(S(1)+p)*sqrt(a+b*x))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when S(2) c d-b e=S(0)

# Integrands of the form (b d+S(2) c d x)**m (a+b x+c x**S(2))**p

# p>S(0)
[(b*d+S(2)*c*d*x)**S(4)*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(40)*(b**S(2)-S(4)*a*c)*d**S(4)*(b+S(2)*c*x)**S(5)/c**S(2)+S(1)/S(56)*d**S(4)*(b+S(2)*c*x)**S(7)/c**S(2)],
[(b*d+S(2)*c*d*x)**S(3)*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(32)*(b**S(2)-S(4)*a*c)*d**S(3)*(b+S(2)*c*x)**S(4)/c**S(2)+S(1)/S(48)*d**S(3)*(b+S(2)*c*x)**S(6)/c**S(2)],
[(b*d+S(2)*c*d*x)**S(2)*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(24)*(b**S(2)-S(4)*a*c)*d**S(2)*(b+S(2)*c*x)**S(3)/c**S(2)+S(1)/S(40)*d**S(2)*(b+S(2)*c*x)**S(5)/c**S(2)],

# p<S(0)
[(b*d+S(2)*c*d*x)**S(8)/(a+b*x+c*x**S(2)),x,S(6),S(2)*(b**S(2)-S(4)*a*c)**S(3)*d**S(8)*(b+S(2)*c*x)+S(2)/S(3)*(b**S(2)-S(4)*a*c)**S(2)*d**S(8)*(b+S(2)*c*x)**S(3)+S(2)/S(5)*(b**S(2)-S(4)*a*c)*d**S(8)*(b+S(2)*c*x)**S(5)+S(2)/S(7)*d**S(8)*(b+S(2)*c*x)**S(7)-S(2)*(b**S(2)-S(4)*a*c)**(S(7)/S(2))*d**S(8)*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))],
[(b*d+S(2)*c*d*x)**S(7)/(a+b*x+c*x**S(2)),x,S(5),(b**S(2)-S(4)*a*c)**S(2)*d**S(7)*(b+S(2)*c*x)**S(2)+S(1)/S(2)*(b**S(2)-S(4)*a*c)*d**S(7)*(b+S(2)*c*x)**S(4)+S(1)/S(3)*d**S(7)*(b+S(2)*c*x)**S(6)+(b**S(2)-S(4)*a*c)**S(3)*d**S(7)*log(a+b*x+c*x**S(2))],
[S(1)/((b*d+S(2)*c*d*x)**S(4)*(a+b*x+c*x**S(2))**S(3)),x,S(6),S(140)/S(3)*c**S(2)/((b**S(2)-S(4)*a*c)**S(3)*d**S(4)*(b+S(2)*c*x)**S(3))+S(140)*c**S(2)/((b**S(2)-S(4)*a*c)**S(4)*d**S(4)*(b+S(2)*c*x))+(-S(1)/S(2))/((b**S(2)-S(4)*a*c)*d**S(4)*(b+S(2)*c*x)**S(3)*(a+b*x+c*x**S(2))**S(2))+S(7)*c/((b**S(2)-S(4)*a*c)**S(2)*d**S(4)*(b+S(2)*c*x)**S(3)*(a+b*x+c*x**S(2)))-S(140)*c**S(2)*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/((b**S(2)-S(4)*a*c)**(S(9)/S(2))*d**S(4))],

# Integrands of the form (b d+S(2) c d x)**m (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[(b*d+S(2)*c*d*x)**S(4)*(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),-S(1)/S(64)*(b**S(2)-S(4)*a*c)**S(3)*d**S(4)*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/c**(S(3)/S(2))-S(1)/S(32)*(b**S(2)-S(4)*a*c)**S(2)*d**S(4)*(b+S(2)*c*x)*sqrt(a+b*x+c*x**S(2))/c-S(1)/S(48)*(b**S(2)-S(4)*a*c)*d**S(4)*(b+S(2)*c*x)**S(3)*sqrt(a+b*x+c*x**S(2))/c+S(1)/S(12)*d**S(4)*(b+S(2)*c*x)**S(5)*sqrt(a+b*x+c*x**S(2))/c],
[(b*d+S(2)*c*d*x)**S(3)*(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(3),S(4)/S(15)*(b**S(2)-S(4)*a*c)*d**S(3)*(a+b*x+c*x**S(2))**(S(3)/S(2))+S(2)/S(5)*d**S(3)*(b+S(2)*c*x)**S(2)*(a+b*x+c*x**S(2))**(S(3)/S(2))],
[(a+b*x+c*x**S(2))**(S(5)/S(2))/(b*d+S(2)*c*d*x)**S(12),x,S(3),S(2)/S(11)*(a+b*x+c*x**S(2))**(S(7)/S(2))/((b**S(2)-S(4)*a*c)*d**S(12)*(b+S(2)*c*x)**S(11))+S(8)/S(99)*(a+b*x+c*x**S(2))**(S(7)/S(2))/((b**S(2)-S(4)*a*c)**S(2)*d**S(12)*(b+S(2)*c*x)**S(9))+S(16)/S(693)*(a+b*x+c*x**S(2))**(S(7)/S(2))/((b**S(2)-S(4)*a*c)**S(3)*d**S(12)*(b+S(2)*c*x)**S(7))],

# p<S(0)
[(b*d+S(2)*c*d*x)**S(4)/(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(4),S(3)/S(8)*(b**S(2)-S(4)*a*c)**S(2)*d**S(4)*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/sqrt(c)+S(3)/S(4)*(b**S(2)-S(4)*a*c)*d**S(4)*(b+S(2)*c*x)*sqrt(a+b*x+c*x**S(2))+S(1)/S(2)*d**S(4)*(b+S(2)*c*x)**S(3)*sqrt(a+b*x+c*x**S(2))],
[(b*d+S(2)*c*d*x)**S(3)/(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(3),S(4)/S(3)*(b**S(2)-S(4)*a*c)*d**S(3)*sqrt(a+b*x+c*x**S(2))+S(2)/S(3)*d**S(3)*(b+S(2)*c*x)**S(2)*sqrt(a+b*x+c*x**S(2))],
[S(1)/((a+b*x)*sqrt(S(1)+a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))),x,S(2),-arctanh(sqrt(S(1)+a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))/b],

# Integrands of the form (b d+S(2) c d x)**(m/S(2)) (a+b x+c x**S(2))**p

# p>S(0)
[(b*d+S(2)*c*d*x)**(S(5)/S(2))*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(28)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(7)/S(2))/(c**S(2)*d)+S(1)/S(44)*(b*d+S(2)*c*d*x)**(S(11)/S(2))/(c**S(2)*d**S(3))],
[(b*d+S(2)*c*d*x)**(S(3)/S(2))*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(20)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(5)/S(2))/(c**S(2)*d)+S(1)/S(36)*(b*d+S(2)*c*d*x)**(S(9)/S(2))/(c**S(2)*d**S(3))],
[(b*d+S(2)*c*d*x)**(S(1)/S(2))*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(12)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(3)/S(2))/(c**S(2)*d)+S(1)/S(28)*(b*d+S(2)*c*d*x)**(S(7)/S(2))/(c**S(2)*d**S(3))],

# p<S(0)
[(b*d+S(2)*c*d*x)**(S(11)/S(2))/(a+b*x+c*x**S(2)),x,S(8),S(4)/S(5)*(b**S(2)-S(4)*a*c)*d**S(3)*(b*d+S(2)*c*d*x)**(S(5)/S(2))+S(4)/S(9)*d*(b*d+S(2)*c*d*x)**(S(9)/S(2))-S(2)*(b**S(2)-S(4)*a*c)**(S(9)/S(4))*d**(S(11)/S(2))*arctan(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))-S(2)*(b**S(2)-S(4)*a*c)**(S(9)/S(4))*d**(S(11)/S(2))*arctanh(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))+S(4)*(b**S(2)-S(4)*a*c)**S(2)*d**S(5)*sqrt(b*d+S(2)*c*d*x)],
[(b*d+S(2)*c*d*x)**(S(9)/S(2))/(a+b*x+c*x**S(2)),x,S(7),S(4)/S(3)*(b**S(2)-S(4)*a*c)*d**S(3)*(b*d+S(2)*c*d*x)**(S(3)/S(2))+S(4)/S(7)*d*(b*d+S(2)*c*d*x)**(S(7)/S(2))+S(2)*(b**S(2)-S(4)*a*c)**(S(7)/S(4))*d**(S(9)/S(2))*arctan(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))-S(2)*(b**S(2)-S(4)*a*c)**(S(7)/S(4))*d**(S(9)/S(2))*arctanh(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))],
[(b*d+S(2)*c*d*x)**(S(7)/S(2))/(a+b*x+c*x**S(2)),x,S(7),S(4)/S(5)*d*(b*d+S(2)*c*d*x)**(S(5)/S(2))-S(2)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*d**(S(7)/S(2))*arctan(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))-S(2)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*d**(S(7)/S(2))*arctanh(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)))+S(4)*(b**S(2)-S(4)*a*c)*d**S(3)*sqrt(b*d+S(2)*c*d*x)],

# Integrands of the form (b d+S(2) c d x)**(m/S(2)) (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[(b*d+S(2)*c*d*x)**(S(7)/S(2))*(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(6),-S(2)/S(77)*(b**S(2)-S(4)*a*c)*d*(b*d+S(2)*c*d*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2))/c+S(1)/S(11)*(b*d+S(2)*c*d*x)**(S(9)/S(2))*sqrt(a+b*x+c*x**S(2))/(c*d)-S(10)/S(231)*(b**S(2)-S(4)*a*c)**S(2)*d**S(3)*sqrt(b*d+S(2)*c*d*x)*sqrt(a+b*x+c*x**S(2))/c-S(5)/S(231)*(b**S(2)-S(4)*a*c)**(S(13)/S(4))*d**(S(7)/S(2))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(2)*sqrt(a+b*x+c*x**S(2)))],
[(b*d+S(2)*c*d*x)**(S(3)/S(2))*(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),S(1)/S(7)*(b*d+S(2)*c*d*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2))/(c*d)-S(2)/S(21)*(b**S(2)-S(4)*a*c)*d*sqrt(b*d+S(2)*c*d*x)*sqrt(a+b*x+c*x**S(2))/c-S(1)/S(21)*(b**S(2)-S(4)*a*c)**(S(9)/S(4))*d**(S(3)/S(2))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(2)*sqrt(a+b*x+c*x**S(2)))],
[(a+b*x+c*x**S(2))**(S(1)/S(2))/(b*d+S(2)*c*d*x)**(S(1)/S(2)),x,S(4),S(1)/S(3)*sqrt(b*d+S(2)*c*d*x)*sqrt(a+b*x+c*x**S(2))/(c*d)-S(1)/S(3)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(2)*sqrt(d)*sqrt(a+b*x+c*x**S(2)))],

# p<S(0)
[(b*d+S(2)*c*d*x)**(S(7)/S(2))/(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(5),S(4)/S(7)*d*(b*d+S(2)*c*d*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2))+S(20)/S(21)*(b**S(2)-S(4)*a*c)*d**S(3)*sqrt(b*d+S(2)*c*d*x)*sqrt(a+b*x+c*x**S(2))+S(10)/S(21)*(b**S(2)-S(4)*a*c)**(S(9)/S(4))*d**(S(7)/S(2))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c*sqrt(a+b*x+c*x**S(2)))],
[(b*d+S(2)*c*d*x)**(S(3)/S(2))/(a+b*x+c*x**S(2))**(S(1)/S(2)),x,S(4),S(4)/S(3)*d*sqrt(b*d+S(2)*c*d*x)*sqrt(a+b*x+c*x**S(2))+S(2)/S(3)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*d**(S(3)/S(2))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c*sqrt(a+b*x+c*x**S(2)))],
[S(1)/((b*d+S(2)*c*d*x)**(S(1)/S(2))*(a+b*x+c*x**S(2))**(S(1)/S(2))),x,S(3),S(2)*(b**S(2)-S(4)*a*c)**(S(1)/S(4))*EllipticF(sqrt(b*d+S(2)*c*d*x)/((b**S(2)-S(4)*a*c)**(S(1)/S(4))*sqrt(d)),I)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c*sqrt(d)*sqrt(a+b*x+c*x**S(2)))],

# Integrands of the form (b d+S(2) c d x)**(m/S(3)) (a+b x+c x**S(2))**(p/S(3))
[(a+b*x+c*x**S(2))**(S(4)/S(3))/(b*d+S(2)*c*d*x)**(S(11)/S(3)),x,S(15),S(3)/S(4)*(a+b*x+c*x**S(2))**(S(7)/S(3))/((b**S(2)-S(4)*a*c)*d*(b*d+S(2)*c*d*x)**(S(8)/S(3)))-S(9)/S(4)*(a+b*x+c*x**S(2))**(S(7)/S(3))/((b**S(2)-S(4)*a*c)**S(2)*d**S(3)*(b*d+S(2)*c*d*x)**(S(2)/S(3)))-S(3)/S(16)*(b*d+S(2)*c*d*x)**(S(4)/S(3))*(S(4)*a-b**S(2)/c+(b+S(2)*c*x)**S(2)/c)**(S(1)/S(3))/(S(2)**(S(2)/S(3))*c**S(2)*(b**S(2)-S(4)*a*c)*d**S(5))+S(9)/S(64)*(b*d+S(2)*c*d*x)**(S(4)/S(3))*(S(4)*a-b**S(2)/c+(b+S(2)*c*x)**S(2)/c)**(S(4)/S(3))/(S(2)**(S(2)/S(3))*c*(b**S(2)-S(4)*a*c)**S(2)*d**S(5))-S(1)/S(16)*log((-S(2)**(S(1)/S(3))*(d*(b+S(2)*c*x))**(S(2)/S(3))+S(2)*c**(S(1)/S(3))*d**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(7)/S(3))*d**(S(11)/S(3)))+S(1)/S(32)*log((b*d*(d*(b+S(2)*c*x))**(S(1)/S(3))+S(2)*c*d*x*(d*(b+S(2)*c*x))**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*d**(S(2)/S(3))*(d*(b+S(2)*c*x))**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*d**(S(4)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(7)/S(3))*d**(S(11)/S(3)))-S(1)/S(16)*arctan((c**(S(1)/S(3))*d**(S(2)/S(3))+S(2)*(b*d+S(2)*c*d*x)**(S(2)/S(3))/(S(4)*a-b**S(2)/c+(b+S(2)*c*x)**S(2)/c)**(S(1)/S(3)))/(c**(S(1)/S(3))*d**(S(2)/S(3))*sqrt(S(3))))*sqrt(S(3))/(S(2)**(S(2)/S(3))*c**(S(7)/S(3))*d**(S(11)/S(3)))],
[(a+b*x+c*x**S(2))**(S(4)/S(3))/(b*d+S(2)*c*d*x)**(S(17)/S(3)),x,S(1),S(3)/S(7)*(a+b*x+c*x**S(2))**(S(7)/S(3))/((b**S(2)-S(4)*a*c)*d*(b*d+S(2)*c*d*x)**(S(14)/S(3)))],

# Integrands of the form (b d+S(2) c d x)**m (a+b x+c x**S(2))**p with m symbolic
[(b*d+S(2)*c*d*x)**m*(a+b*x+c*x**S(2))**S(3),x,S(2),-S(1)/S(128)*(b**S(2)-S(4)*a*c)**S(3)*(b*d+S(2)*c*d*x)**(S(1)+m)/(c**S(4)*d*(S(1)+m))+S(3)/S(128)*(b**S(2)-S(4)*a*c)**S(2)*(b*d+S(2)*c*d*x)**(S(3)+m)/(c**S(4)*d**S(3)*(S(3)+m))-S(3)/S(128)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(5)+m)/(c**S(4)*d**S(5)*(S(5)+m))+S(1)/S(128)*(b*d+S(2)*c*d*x)**(S(7)+m)/(c**S(4)*d**S(7)*(S(7)+m))],
[(b*d+S(2)*c*d*x)**m*(a+b*x+c*x**S(2))**S(2),x,S(2),S(1)/S(32)*(b**S(2)-S(4)*a*c)**S(2)*(b*d+S(2)*c*d*x)**(S(1)+m)/(c**S(3)*d*(S(1)+m))-S(1)/S(16)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(3)+m)/(c**S(3)*d**S(3)*(S(3)+m))+S(1)/S(32)*(b*d+S(2)*c*d*x)**(S(5)+m)/(c**S(3)*d**S(5)*(S(5)+m))],
[(b*d+S(2)*c*d*x)**m*(a+b*x+c*x**S(2)),x,S(2),-S(1)/S(8)*(b**S(2)-S(4)*a*c)*(b*d+S(2)*c*d*x)**(S(1)+m)/(c**S(2)*d*(S(1)+m))+S(1)/S(8)*(b*d+S(2)*c*d*x)**(S(3)+m)/(c**S(2)*d**S(3)*(S(3)+m))],

# Integrands of the form (b d+S(2) c d x)**m (a+b x+c x**S(2))**p with p symbolic
[(b*d+S(2)*c*d*x)**m*(a+b*x+c*x**S(2))**p,x,S(3),-S(2)*(b*d+S(2)*c*d*x)**(S(1)+m)*(S(1)/S(4)*(S(4)*a-b**S(2)/c)+S(1)/S(4)*(b+S(2)*c*x)**S(2)/c)**(S(1)+p)*hypergeom([S(1),S(1)/S(2)*(S(3)+m+S(2)*p)],[S(1)/S(2)*(S(3)+m)],(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/((b**S(2)-S(4)*a*c)*d*(S(1)+m)),S(1)/S(2)*(b*d+S(2)*c*d*x)**(S(1)+m)*(S(1)/S(4)*(S(4)*a-b**S(2)/c)+S(1)/S(4)*(b+S(2)*c*x)**S(2)/c)**p*hypergeom([S(1)/S(2)*(S(1)+m),-p],[S(1)/S(2)*(S(3)+m)],(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/(c*d*(S(1)+m)*(S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**p)],
[(b*d+S(2)*c*d*x)**S(5)*(a+b*x+c*x**S(2))**p,x,S(4),S(2)*(b**S(2)-S(4)*a*c)**S(2)*d**S(5)*(a+b*x+c*x**S(2))**(S(1)+p)/((S(3)+p)*(S(2)+S(3)*p+p**S(2)))+S(2)*(b**S(2)-S(4)*a*c)*d**S(5)*(b+S(2)*c*x)**S(2)*(a+b*x+c*x**S(2))**(S(1)+p)/(S(6)+S(5)*p+p**S(2))+d**S(5)*(b+S(2)*c*x)**S(4)*(a+b*x+c*x**S(2))**(S(1)+p)/(S(3)+p)],
[(b*d+S(2)*c*d*x)**S(4)*(a+b*x+c*x**S(2))**p,x,S(3),-S(2)/S(5)*d**S(4)*(b+S(2)*c*x)**S(5)*(S(1)/S(4)*(S(4)*a-b**S(2)/c)+S(1)/S(4)*(b+S(2)*c*x)**S(2)/c)**(S(1)+p)*hypergeom([S(1),S(7)/S(2)+p],[S(7)/S(2)],(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/(b**S(2)-S(4)*a*c),S(1)/S(10)*d**S(4)*(b+S(2)*c*x)**S(5)*(S(1)/S(4)*(S(4)*a-b**S(2)/c)+S(1)/S(4)*(b+S(2)*c*x)**S(2)/c)**p*hypergeom([S(5)/S(2),-p],[S(7)/S(2)],(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/(c*(S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**p)],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when b**S(2)-S(4) a c=S(0)

# Integrands of the form (d+e x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p

# p>S(0)
[(d+e*x)**S(4)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(1)/S(5)*(b*d-a*e)**S(2)*(d+e*x)**S(5)/e**S(3)-S(1)/S(3)*b*(b*d-a*e)*(d+e*x)**S(6)/e**S(3)+S(1)/S(7)*b**S(2)*(d+e*x)**S(7)/e**S(3)],
[(d+e*x)**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(1)/S(4)*(b*d-a*e)**S(2)*(d+e*x)**S(4)/e**S(3)-S(2)/S(5)*b*(b*d-a*e)*(d+e*x)**S(5)/e**S(3)+S(1)/S(6)*b**S(2)*(d+e*x)**S(6)/e**S(3)],
[(d+e*x)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(1)/S(3)*(b*d-a*e)**S(2)*(a+b*x)**S(3)/b**S(3)+S(1)/S(2)*e*(b*d-a*e)*(a+b*x)**S(4)/b**S(3)+S(1)/S(5)*e**S(2)*(a+b*x)**S(5)/b**S(3)],

# p<S(0)
[(d+e*x)**S(5)/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(10)*e**S(2)*(b*d-a*e)**S(3)*x/b**S(5)-(b*d-a*e)**S(5)/(b**S(6)*(a+b*x))+S(5)*e**S(3)*(b*d-a*e)**S(2)*(a+b*x)**S(2)/b**S(6)+S(5)/S(3)*e**S(4)*(b*d-a*e)*(a+b*x)**S(3)/b**S(6)+S(1)/S(4)*e**S(5)*(a+b*x)**S(4)/b**S(6)+S(5)*e*(b*d-a*e)**S(4)*log(a+b*x)/b**S(6)],
[(d+e*x)/(S(9)+S(12)*x+S(4)*x**S(2))**S(2),x,S(3),S(1)/S(12)*(-S(2)*d+S(3)*e)/(S(3)+S(2)*x)**S(3)-S(1)/S(8)*e/(S(3)+S(2)*x)**S(2)],
[(d+e*x)/(S(9)+S(12)*x+S(4)*x**S(2))**S(3),x,S(3),S(1)/S(20)*(-S(2)*d+S(3)*e)/(S(3)+S(2)*x)**S(5)-S(1)/S(16)*e/(S(3)+S(2)*x)**S(4)],

# Integrands of the form (d+e x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**S(4)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(1)/S(6)*(d+e*x)**S(5)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(1)/S(30)*(b*d-a*e)*(d+e*x)**S(5)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],
[(d+e*x)**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(1)/S(5)*(d+e*x)**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(1)/S(20)*(b*d-a*e)*(d+e*x)**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],
[(d+e*x)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(1)/S(4)*(d+e*x)**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(1)/S(12)*(b*d-a*e)*(d+e*x)**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],
[(d+e*x)*(S(9)+S(12)*x+S(4)*x**S(2))**(S(1)/S(2)),x,S(2),S(1)/S(12)*e*(S(9)+S(12)*x+S(4)*x**S(2))**(S(3)/S(2))+S(1)/S(8)*(S(2)*d-S(3)*e)*(S(3)+S(2)*x)*sqrt(S(9)+S(12)*x+S(4)*x**S(2))],

# p<S(0)
[(d+e*x)**S(4)/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(6),S(1)/S(2)*(b*d-a*e)**S(2)*(a+b*x)*(d+e*x)**S(2)/(b**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(1)/S(3)*(b*d-a*e)*(a+b*x)*(d+e*x)**S(3)/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(1)/S(4)*(a+b*x)*(d+e*x)**S(4)/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+(b*d-a*e)**S(4)*(a+b*x)*log(a+b*x)/(b**S(5)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+e*(b*d-a*e)**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/b**S(5)],
[(d+e*x)**S(3)/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(5),S(1)/S(2)*(b*d-a*e)*(a+b*x)*(d+e*x)**S(2)/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(1)/S(3)*(a+b*x)*(d+e*x)**S(3)/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+(b*d-a*e)**S(3)*(a+b*x)*log(a+b*x)/(b**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+e*(b*d-a*e)**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/b**S(4)],
[(d+e*x)**S(2)/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(4),S(1)/S(2)*(a+b*x)*(d+e*x)**S(2)/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+(b*d-a*e)**S(2)*(a+b*x)*log(a+b*x)/(b**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+e*(b*d-a*e)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/b**S(3)],

# Integrands of the form (d+e x)**(m/S(2)) (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p

# p>S(0)
[(d+e*x)**(S(7)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(2)/S(9)*(b*d-a*e)**S(2)*(d+e*x)**(S(9)/S(2))/e**S(3)-S(4)/S(11)*b*(b*d-a*e)*(d+e*x)**(S(11)/S(2))/e**S(3)+S(2)/S(13)*b**S(2)*(d+e*x)**(S(13)/S(2))/e**S(3)],
[(d+e*x)**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(2)/S(7)*(b*d-a*e)**S(2)*(d+e*x)**(S(7)/S(2))/e**S(3)-S(4)/S(9)*b*(b*d-a*e)*(d+e*x)**(S(9)/S(2))/e**S(3)+S(2)/S(11)*b**S(2)*(d+e*x)**(S(11)/S(2))/e**S(3)],
[(d+e*x)**(S(3)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),S(2)/S(5)*(b*d-a*e)**S(2)*(d+e*x)**(S(5)/S(2))/e**S(3)-S(4)/S(7)*b*(b*d-a*e)*(d+e*x)**(S(7)/S(2))/e**S(3)+S(2)/S(9)*b**S(2)*(d+e*x)**(S(9)/S(2))/e**S(3)],

# p<S(0)
[(d+e*x)**(S(9)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(8),S(3)*e*(b*d-a*e)**S(2)*(d+e*x)**(S(3)/S(2))/b**S(4)+S(9)/S(5)*e*(b*d-a*e)*(d+e*x)**(S(5)/S(2))/b**S(3)+S(9)/S(7)*e*(d+e*x)**(S(7)/S(2))/b**S(2)-(d+e*x)**(S(9)/S(2))/(b*(a+b*x))-S(9)*e*(b*d-a*e)**(S(7)/S(2))*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/b**(S(11)/S(2))+S(9)*e*(b*d-a*e)**S(3)*sqrt(d+e*x)/b**S(5)],
[(d+e*x)**(S(7)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(7),S(7)/S(3)*e*(b*d-a*e)*(d+e*x)**(S(3)/S(2))/b**S(3)+S(7)/S(5)*e*(d+e*x)**(S(5)/S(2))/b**S(2)-(d+e*x)**(S(7)/S(2))/(b*(a+b*x))-S(7)*e*(b*d-a*e)**(S(5)/S(2))*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/b**(S(9)/S(2))+S(7)*e*(b*d-a*e)**S(2)*sqrt(d+e*x)/b**S(4)],
[(d+e*x)**(S(5)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(6),S(5)/S(3)*e*(d+e*x)**(S(3)/S(2))/b**S(2)-(d+e*x)**(S(5)/S(2))/(b*(a+b*x))-S(5)*e*(b*d-a*e)**(S(3)/S(2))*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/b**(S(7)/S(2))+S(5)*e*(b*d-a*e)*sqrt(d+e*x)/b**S(3)],

# Integrands of the form (d+e x)**(m/S(2)) (a**S(2)+S(2) a b x+b**S(2) x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**(S(5)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(2)/S(9)*(d+e*x)**(S(7)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(4)/S(63)*(b*d-a*e)*(d+e*x)**(S(7)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],
[(d+e*x)**(S(3)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(2)/S(7)*(d+e*x)**(S(5)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(4)/S(35)*(b*d-a*e)*(d+e*x)**(S(5)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],
[(d+e*x)**(S(1)/S(2))*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(2),S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/e-S(4)/S(15)*(b*d-a*e)*(d+e*x)**(S(3)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))/(e**S(2)*(a+b*x))],

# p<S(0)
[(d+e*x)**(S(7)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(7),S(2)/S(3)*(b*d-a*e)**S(2)*(a+b*x)*(d+e*x)**(S(3)/S(2))/(b**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)/S(5)*(b*d-a*e)*(a+b*x)*(d+e*x)**(S(5)/S(2))/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)/S(7)*(a+b*x)*(d+e*x)**(S(7)/S(2))/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(2)*(b*d-a*e)**(S(7)/S(2))*(a+b*x)*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/(b**(S(9)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)*(b*d-a*e)**S(3)*(a+b*x)*sqrt(d+e*x)/(b**S(4)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],
[(d+e*x)**(S(5)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(6),S(2)/S(3)*(b*d-a*e)*(a+b*x)*(d+e*x)**(S(3)/S(2))/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)/S(5)*(a+b*x)*(d+e*x)**(S(5)/S(2))/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(2)*(b*d-a*e)**(S(5)/S(2))*(a+b*x)*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/(b**(S(7)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)*(b*d-a*e)**S(2)*(a+b*x)*sqrt(d+e*x)/(b**S(3)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],
[(d+e*x)**(S(3)/S(2))/(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)/S(2)),x,S(5),S(2)/S(3)*(a+b*x)*(d+e*x)**(S(3)/S(2))/(b*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))-S(2)*(b*d-a*e)**(S(3)/S(2))*(a+b*x)*arctanh(sqrt(b)*sqrt(d+e*x)/sqrt(b*d-a*e))/(b**(S(5)/S(2))*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))+S(2)*(b*d-a*e)*(a+b*x)*sqrt(d+e*x)/(b**S(2)*sqrt(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)))],

# Integrands of the form (d+e x)**m (c d**S(2)+S(2) c d e x+c e**S(2) x**S(2))**p

# p>S(0)
[(d+e*x)**m*(c*d**S(2)+S(2)*c*d*e*x+c*e**S(2)*x**S(2)),x,S(3),c*(d+e*x)**(S(3)+m)/(e*(S(3)+m))],
[(d+e*x)**S(2)*(c*d**S(2)+S(2)*c*d*e*x+c*e**S(2)*x**S(2)),x,S(3),S(1)/S(5)*c*(d+e*x)**S(5)/e],

# p<S(0)
[(d+e*x)**m/(c*d**S(2)+S(2)*c*d*e*x+c*e**S(2)*x**S(2)),x,S(3),-(d+e*x)**(-S(1)+m)/(c*e*(S(1)-m))],
[(d+e*x)**S(5)/(c*d**S(2)+S(2)*c*d*e*x+c*e**S(2)*x**S(2)),x,S(3),S(1)/S(4)*(d+e*x)**S(4)/(c*e)],
[(d+e*x)**S(4)/(c*d**S(2)+S(2)*c*d*e*x+c*e**S(2)*x**S(2)),x,S(3),S(1)/S(3)*(d+e*x)**S(3)/(c*e)],

# Integrands of the form (d+e x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p with m symbolic
[(d+e*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(3),x,S(3),(b*d-a*e)**S(6)*(d+e*x)**(S(1)+m)/(e**S(7)*(S(1)+m))-S(6)*b*(b*d-a*e)**S(5)*(d+e*x)**(S(2)+m)/(e**S(7)*(S(2)+m))+S(15)*b**S(2)*(b*d-a*e)**S(4)*(d+e*x)**(S(3)+m)/(e**S(7)*(S(3)+m))-S(20)*b**S(3)*(b*d-a*e)**S(3)*(d+e*x)**(S(4)+m)/(e**S(7)*(S(4)+m))+S(15)*b**S(4)*(b*d-a*e)**S(2)*(d+e*x)**(S(5)+m)/(e**S(7)*(S(5)+m))-S(6)*b**S(5)*(b*d-a*e)*(d+e*x)**(S(6)+m)/(e**S(7)*(S(6)+m))+b**S(6)*(d+e*x)**(S(7)+m)/(e**S(7)*(S(7)+m))],
[(d+e*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**S(2),x,S(3),(b*d-a*e)**S(4)*(d+e*x)**(S(1)+m)/(e**S(5)*(S(1)+m))-S(4)*b*(b*d-a*e)**S(3)*(d+e*x)**(S(2)+m)/(e**S(5)*(S(2)+m))+S(6)*b**S(2)*(b*d-a*e)**S(2)*(d+e*x)**(S(3)+m)/(e**S(5)*(S(3)+m))-S(4)*b**S(3)*(b*d-a*e)*(d+e*x)**(S(4)+m)/(e**S(5)*(S(4)+m))+b**S(4)*(d+e*x)**(S(5)+m)/(e**S(5)*(S(5)+m))],
[(d+e*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2)),x,S(3),(b*d-a*e)**S(2)*(d+e*x)**(S(1)+m)/(e**S(3)*(S(1)+m))-S(2)*b*(b*d-a*e)*(d+e*x)**(S(2)+m)/(e**S(3)*(S(2)+m))+b**S(2)*(d+e*x)**(S(3)+m)/(e**S(3)*(S(3)+m))],

# Integrands of the form (d+e x)**m (a**S(2)+S(2) a b x+b**S(2) x**S(2))**p with p symbolic
[(d+e*x)**m*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p,x,S(3),(d+e*x)**(S(1)+m)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p*hypergeom([S(1)+m,-S(2)*p],[S(2)+m],b*(d+e*x)/(b*d-a*e))/(e*(S(1)+m)*(-e*(a+b*x)/(b*d-a*e))**(S(2)*p))],
[(d+e*x)**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p,x,S(4),S(3)*(b*d-a*e)**S(3)*(a+b*x)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b**S(4)*(S(2)+p)*(S(3)+S(8)*p+S(4)*p**S(2)))+S(3)/S(2)*(b*d-a*e)*(a+b*x)*(d+e*x)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b**S(2)*(S(6)+S(7)*p+S(2)*p**S(2)))+S(1)/S(2)*(a+b*x)*(d+e*x)**S(3)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b*(S(2)+p))+S(3)/S(2)*e*(b*d-a*e)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)+p)/(b**S(4)*(S(2)+p)*(S(3)+S(5)*p+S(2)*p**S(2)))],
[(d+e*x)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p,x,S(3),S(2)*(b*d-a*e)**S(2)*(a+b*x)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b**S(3)*(S(3)+S(8)*p+S(4)*p**S(2)))+(a+b*x)*(d+e*x)**S(2)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**p/(b*(S(3)+S(2)*p))+e*(b*d-a*e)*(a**S(2)+S(2)*a*b*x+b**S(2)*x**S(2))**(S(1)+p)/(b**S(3)*(S(1)+p)*(S(3)+S(2)*p))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p when c d**S(2)-b d e+a e**S(2)=S(0)

# Integrands of the form (d+e x)**m (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**p

# p>S(0)
[(d+e*x)**S(4)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),S(1)/S(6)*(a-c*d**S(2)/e**S(2))*(d+e*x)**S(6)+S(1)/S(7)*c*d*(d+e*x)**S(7)/e**S(2)],
[(d+e*x)**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),S(1)/S(5)*(a-c*d**S(2)/e**S(2))*(d+e*x)**S(5)+S(1)/S(6)*c*d*(d+e*x)**S(6)/e**S(2)],
[(d+e*x)**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),S(1)/S(4)*(a-c*d**S(2)/e**S(2))*(d+e*x)**S(4)+S(1)/S(5)*c*d*(d+e*x)**S(5)/e**S(2)],

# p<S(0)
[(d+e*x)**S(5)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),e*(c*d**S(2)-a*e**S(2))**S(3)*x/(c**S(4)*d**S(4))+S(1)/S(2)*(c*d**S(2)-a*e**S(2))**S(2)*(d+e*x)**S(2)/(c**S(3)*d**S(3))+S(1)/S(3)*(c*d**S(2)-a*e**S(2))*(d+e*x)**S(3)/(c**S(2)*d**S(2))+S(1)/S(4)*(d+e*x)**S(4)/(c*d)+(c*d**S(2)-a*e**S(2))**S(4)*log(a*e+c*d*x)/(c**S(5)*d**S(5))],
[(d+e*x)**S(4)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),e*(c*d**S(2)-a*e**S(2))**S(2)*x/(c**S(3)*d**S(3))+S(1)/S(2)*(c*d**S(2)-a*e**S(2))*(d+e*x)**S(2)/(c**S(2)*d**S(2))+S(1)/S(3)*(d+e*x)**S(3)/(c*d)+(c*d**S(2)-a*e**S(2))**S(3)*log(a*e+c*d*x)/(c**S(4)*d**S(4))],
[(d+e*x)**S(3)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),e*(c*d**S(2)-a*e**S(2))*x/(c**S(2)*d**S(2))+S(1)/S(2)*(d+e*x)**S(2)/(c*d)+(c*d**S(2)-a*e**S(2))**S(2)*log(a*e+c*d*x)/(c**S(3)*d**S(3))],

# Integrands of the form (d+e x)**m (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**S(4)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(7),S(7)/S(64)*(c*d**S(2)-a*e**S(2))**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(4)*d**S(4))+S(21)/S(160)*(c*d**S(2)-a*e**S(2))**S(2)*(d+e*x)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(3)*d**S(3))+S(3)/S(20)*(c*d**S(2)-a*e**S(2))*(d+e*x)**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(2)*d**S(2))+S(1)/S(6)*(d+e*x)**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c*d)-S(21)/S(1024)*(c*d**S(2)-a*e**S(2))**S(6)*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(11)/S(2))*d**(S(11)/S(2))*e**(S(3)/S(2)))+S(21)/S(512)*(c*d**S(2)-a*e**S(2))**S(4)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(5)*d**S(5)*e)],
[(d+e*x)**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(6),S(7)/S(48)*(c*d**S(2)-a*e**S(2))**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(3)*d**S(3))+S(7)/S(40)*(c*d**S(2)-a*e**S(2))*(d+e*x)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(2)*d**S(2))+S(1)/S(5)*(d+e*x)**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c*d)-S(7)/S(256)*(c*d**S(2)-a*e**S(2))**S(5)*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(9)/S(2))*d**(S(9)/S(2))*e**(S(3)/S(2)))+S(7)/S(128)*(c*d**S(2)-a*e**S(2))**S(3)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(4)*d**S(4)*e)],
[(d+e*x)**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(5),S(5)/S(24)*(c*d**S(2)-a*e**S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(2)*d**S(2))+S(1)/S(4)*(d+e*x)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c*d)-S(5)/S(128)*(c*d**S(2)-a*e**S(2))**S(4)*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(7)/S(2))*d**(S(7)/S(2))*e**(S(3)/S(2)))+S(5)/S(64)*(c*d**S(2)-a*e**S(2))**S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(3)*d**S(3)*e)],

# p<S(0)
[(d+e*x)**S(3)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(5),S(5)/S(16)*(c*d**S(2)-a*e**S(2))**S(3)*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(7)/S(2))*d**(S(7)/S(2))*sqrt(e))+S(5)/S(8)*(c*d**S(2)-a*e**S(2))**S(2)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(3)*d**S(3))+S(5)/S(12)*(c*d**S(2)-a*e**S(2))*(d+e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(2)*d**S(2))+S(1)/S(3)*(d+e*x)**S(2)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c*d)],
[(d+e*x)**S(2)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(4),S(3)/S(8)*(c*d**S(2)-a*e**S(2))**S(2)*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(5)/S(2))*d**(S(5)/S(2))*sqrt(e))+S(3)/S(4)*(c*d**S(2)-a*e**S(2))*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(2)*d**S(2))+S(1)/S(2)*(d+e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c*d)],
[(d+e*x)/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(3),S(1)/S(2)*(c*d**S(2)-a*e**S(2))*arctanh(S(1)/S(2)*(c*d**S(2)+a*e**S(2)+S(2)*c*d*e*x)/(sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))))/(c**(S(3)/S(2))*d**(S(3)/S(2))*sqrt(e))+sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c*d)],

# Integrands of the form (d+e x)**(m/S(2)) (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**p

# p>S(0)
[(d+e*x)**(S(3)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),S(2)/S(7)*(a-c*d**S(2)/e**S(2))*(d+e*x)**(S(7)/S(2))+S(2)/S(9)*c*d*(d+e*x)**(S(9)/S(2))/e**S(2)],
[(d+e*x)**(S(1)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(3),S(2)/S(5)*(a-c*d**S(2)/e**S(2))*(d+e*x)**(S(5)/S(2))+S(2)/S(7)*c*d*(d+e*x)**(S(7)/S(2))/e**S(2)],
[(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(d+e*x)**(S(1)/S(2)),x,S(3),S(2)/S(3)*(a-c*d**S(2)/e**S(2))*(d+e*x)**(S(3)/S(2))+S(2)/S(5)*c*d*(d+e*x)**(S(5)/S(2))/e**S(2)],

# p<S(0)
[(d+e*x)**(S(9)/S(2))/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(7),S(2)/S(3)*(c*d**S(2)-a*e**S(2))**S(2)*(d+e*x)**(S(3)/S(2))/(c**S(3)*d**S(3))+S(2)/S(5)*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(5)/S(2))/(c**S(2)*d**S(2))+S(2)/S(7)*(d+e*x)**(S(7)/S(2))/(c*d)-S(2)*(c*d**S(2)-a*e**S(2))**(S(7)/S(2))*arctanh(sqrt(c)*sqrt(d)*sqrt(d+e*x)/sqrt(c*d**S(2)-a*e**S(2)))/(c**(S(9)/S(2))*d**(S(9)/S(2)))+S(2)*(c*d**S(2)-a*e**S(2))**S(3)*sqrt(d+e*x)/(c**S(4)*d**S(4))],
[(d+e*x)**(S(7)/S(2))/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(6),S(2)/S(3)*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(3)/S(2))/(c**S(2)*d**S(2))+S(2)/S(5)*(d+e*x)**(S(5)/S(2))/(c*d)-S(2)*(c*d**S(2)-a*e**S(2))**(S(5)/S(2))*arctanh(sqrt(c)*sqrt(d)*sqrt(d+e*x)/sqrt(c*d**S(2)-a*e**S(2)))/(c**(S(7)/S(2))*d**(S(7)/S(2)))+S(2)*(c*d**S(2)-a*e**S(2))**S(2)*sqrt(d+e*x)/(c**S(3)*d**S(3))],
[(d+e*x)**(S(5)/S(2))/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(5),S(2)/S(3)*(d+e*x)**(S(3)/S(2))/(c*d)-S(2)*(c*d**S(2)-a*e**S(2))**(S(3)/S(2))*arctanh(sqrt(c)*sqrt(d)*sqrt(d+e*x)/sqrt(c*d**S(2)-a*e**S(2)))/(c**(S(5)/S(2))*d**(S(5)/S(2)))+S(2)*(c*d**S(2)-a*e**S(2))*sqrt(d+e*x)/(c**S(2)*d**S(2))],

# Integrands of the form (d+e x)**(m/S(2)) (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**(S(7)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(5),S(256)/S(3465)*(c*d**S(2)-a*e**S(2))**S(4)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(5)*d**S(5)*(d+e*x)**(S(3)/S(2)))+S(16)/S(99)*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(3)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(2)*d**S(2))+S(2)/S(11)*(d+e*x)**(S(5)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c*d)+S(128)/S(1155)*(c*d**S(2)-a*e**S(2))**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(4)*d**S(4)*sqrt(d+e*x))+S(32)/S(231)*(c*d**S(2)-a*e**S(2))**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/(c**S(3)*d**S(3))],
[(d+e*x)**(S(5)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(4),S(32)/S(315)*(c*d**S(2)-a*e**S(2))**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(4)*d**S(4)*(d+e*x)**(S(3)/S(2)))+S(2)/S(9)*(d+e*x)**(S(3)/S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c*d)+S(16)/S(105)*(c*d**S(2)-a*e**S(2))**S(2)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))/(c**S(3)*d**S(3)*sqrt(d+e*x))+S(4)/S(21)*(c*d**S(2)-a*e**S(2))*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/(c**S(2)*d**S(2))],

# p<S(0)
[(d+e*x)**(S(7)/S(2))/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(4),S(12)/S(35)*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(3)/S(2))*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(2)*d**S(2))+S(2)/S(7)*(d+e*x)**(S(5)/S(2))*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c*d)+S(32)/S(35)*(c*d**S(2)-a*e**S(2))**S(3)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(4)*d**S(4)*sqrt(d+e*x))+S(16)/S(35)*(c*d**S(2)-a*e**S(2))**S(2)*sqrt(d+e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(3)*d**S(3))],
[(d+e*x)**(S(5)/S(2))/(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**(S(1)/S(2)),x,S(3),S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c*d)+S(16)/S(15)*(c*d**S(2)-a*e**S(2))**S(2)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(3)*d**S(3)*sqrt(d+e*x))+S(8)/S(15)*(c*d**S(2)-a*e**S(2))*sqrt(d+e*x)*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))/(c**S(2)*d**S(2))],

# Integrands of the form (d+e x)**(m/S(3)) (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**(p/S(2))
[(d+e*x)**(S(2)/S(3))/sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)),x,S(4),S(3)/S(2)*(a*e+c*d*x)*(d+e*x)**(S(2)/S(3))/(c*d*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2)))+S(1)/S(4)*S(3)**(S(3)/S(4))*(c*d**S(2)-a*e**S(2))**(S(2)/S(3))*(d+e*x)**(S(2)/S(3))*((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3)))*sqrt(cos(arccos(((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)+sqrt(S(3))))))**S(2))/cos(arccos(((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)+sqrt(S(3))))))*EllipticF(sin(arccos(((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)+sqrt(S(3)))))),sqrt(S(1)/S(4)*(S(2)+sqrt(S(3)))))*sqrt(a*e+c*d*x)*sqrt(((c*d**S(2)-a*e**S(2))**(S(2)/S(3))+c**(S(1)/S(3))*d**(S(1)/S(3))*(c*d**S(2)-a*e**S(2))**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))+c**(S(2)/S(3))*d**(S(2)/S(3))*(d+e*x)**(S(2)/S(3)))/((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(c*d*e*sqrt(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))*sqrt(-c*d**S(2)/e+a*e+c*d*(d+e*x)/e)*sqrt(-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3)))/((c*d**S(2)-a*e**S(2))**(S(1)/S(3))-c**(S(1)/S(3))*d**(S(1)/S(3))*(d+e*x)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))],

# Integrands of the form (a+b x)**m (a c+(b c+a d) x+b d x**S(2))**p

# p>S(0)
[(a+b*x)**m*(a*c+(b*c+a*d)*x+b*d*x**S(2)),x,S(3),(b*c-a*d)*(a+b*x)**(S(2)+m)/(b**S(2)*(S(2)+m))+d*(a+b*x)**(S(3)+m)/(b**S(2)*(S(3)+m))],
[(a+b*x)**S(3)*(a*c+(b*c+a*d)*x+b*d*x**S(2)),x,S(3),S(1)/S(5)*(b*c-a*d)*(a+b*x)**S(5)/b**S(2)+S(1)/S(6)*d*(a+b*x)**S(6)/b**S(2)],
[(a+b*x)**S(2)*(a*c+(b*c+a*d)*x+b*d*x**S(2)),x,S(3),S(1)/S(4)*(b*c-a*d)*(a+b*x)**S(4)/b**S(2)+S(1)/S(5)*d*(a+b*x)**S(5)/b**S(2)],

# p<S(0)
[(a+b*x)**S(6)/(a*c+(b*c+a*d)*x+b*d*x**S(2)),x,S(3),b*(b*c-a*d)**S(4)*x/d**S(5)-S(1)/S(2)*(b*c-a*d)**S(3)*(a+b*x)**S(2)/d**S(4)+S(1)/S(3)*(b*c-a*d)**S(2)*(a+b*x)**S(3)/d**S(3)-S(1)/S(4)*(b*c-a*d)*(a+b*x)**S(4)/d**S(2)+S(1)/S(5)*(a+b*x)**S(5)/d-(b*c-a*d)**S(5)*log(c+d*x)/d**S(6)],
[(a+b*x)**S(5)/(a*c+(b*c+a*d)*x+b*d*x**S(2)),x,S(3),-b*(b*c-a*d)**S(3)*x/d**S(4)+S(1)/S(2)*(b*c-a*d)**S(2)*(a+b*x)**S(2)/d**S(3)-S(1)/S(3)*(b*c-a*d)*(a+b*x)**S(3)/d**S(2)+S(1)/S(4)*(a+b*x)**S(4)/d+(b*c-a*d)**S(4)*log(c+d*x)/d**S(5)],

# Integrands of the form (d+e x)**m (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**p with m symbolic
[(d+e*x)**m*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**S(3),x,S(3),-(c*d**S(2)-a*e**S(2))**S(3)*(d+e*x)**(S(4)+m)/(e**S(4)*(S(4)+m))+S(3)*c*d*(c*d**S(2)-a*e**S(2))**S(2)*(d+e*x)**(S(5)+m)/(e**S(4)*(S(5)+m))-S(3)*c**S(2)*d**S(2)*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(6)+m)/(e**S(4)*(S(6)+m))+c**S(3)*d**S(3)*(d+e*x)**(S(7)+m)/(e**S(4)*(S(7)+m))],
[(d+e*x)**m*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**S(2),x,S(3),(c*d**S(2)-a*e**S(2))**S(2)*(d+e*x)**(S(3)+m)/(e**S(3)*(S(3)+m))-S(2)*c*d*(c*d**S(2)-a*e**S(2))*(d+e*x)**(S(4)+m)/(e**S(3)*(S(4)+m))+c**S(2)*d**S(2)*(d+e*x)**(S(5)+m)/(e**S(3)*(S(5)+m))],

# Integrands of the form (d+e x)**m (a d e+(c d**S(2)+a e**S(2))x+c d e x**S(2))**p with p symbolic
[(d+e*x)**m*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p,x,S(3),-(a*e+c*d*x)*(d+e*x)**(S(1)+m)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p*hypergeom([S(1),S(2)+m+S(2)*p],[S(2)+m+p],c*d*(d+e*x)/(c*d**S(2)-a*e**S(2)))/((c*d**S(2)-a*e**S(2))*(S(1)+m+p)),(a*e+c*d*x)*(d+e*x)**m*(c*d*(d+e*x)/(c*d**S(2)-a*e**S(2)))**(-m-p)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p*hypergeom([-m-p,S(1)+p],[S(2)+p],-e*(a*e+c*d*x)/(c*d**S(2)-a*e**S(2)))/(c*d*(S(1)+p))],
[(d+e*x)**S(3)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p,x,S(3),-(a*e+c*d*x)*(d+e*x)**S(4)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p*hypergeom([S(1),S(5)+S(2)*p],[S(5)+p],c*d*(d+e*x)/(c*d**S(2)-a*e**S(2)))/((c*d**S(2)-a*e**S(2))*(S(4)+p)),(c*d**S(2)-a*e**S(2))**S(3)*(a*e+c*d*x)*(a*d*e+(c*d**S(2)+a*e**S(2))*x+c*d*e*x**S(2))**p*hypergeom([-S(3)-p,S(1)+p],[S(2)+p],-e*(a*e+c*d*x)/(c*d**S(2)-a*e**S(2)))/(c**S(4)*d**S(4)*(S(1)+p)*(c*d*(d+e*x)/(c*d**S(2)-a*e**S(2)))**p)],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p

# p>S(0)
[(d+e*x)**S(4)*(a+b*x+c*x**S(2)),x,S(2),S(1)/S(5)*(c*d**S(2)-b*d*e+a*e**S(2))*(d+e*x)**S(5)/e**S(3)-S(1)/S(6)*(S(2)*c*d-b*e)*(d+e*x)**S(6)/e**S(3)+S(1)/S(7)*c*(d+e*x)**S(7)/e**S(3)],
[(d+e*x)**S(3)*(a+b*x+c*x**S(2)),x,S(2),S(1)/S(4)*(c*d**S(2)-b*d*e+a*e**S(2))*(d+e*x)**S(4)/e**S(3)-S(1)/S(5)*(S(2)*c*d-b*e)*(d+e*x)**S(5)/e**S(3)+S(1)/S(6)*c*(d+e*x)**S(6)/e**S(3)],
[(d+e*x)**S(2)*(a+b*x+c*x**S(2)),x,S(2),a*d**S(2)*x+S(1)/S(2)*d*(b*d+S(2)*a*e)*x**S(2)+S(1)/S(3)*(c*d**S(2)+e*(S(2)*b*d+a*e))*x**S(3)+S(1)/S(4)*e*(S(2)*c*d+b*e)*x**S(4)+S(1)/S(5)*c*e**S(2)*x**S(5)],
[(d+e*x)*(a+b*x+c*x**S(2)),x,S(2),a*d*x+S(1)/S(2)*(b*d+a*e)*x**S(2)+S(1)/S(3)*(c*d+b*e)*x**S(3)+S(1)/S(4)*c*e*x**S(4)],
[a+b*x+c*x**S(2),x,S(1),a*x+S(1)/S(2)*b*x**S(2)+S(1)/S(3)*c*x**S(3)],

# p<S(0)
[(d+e*x)**S(4)/(a+b*x+c*x**S(2)),x,S(7),e**S(2)*(S(6)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(4)*b*d+a*e))*x/c**S(3)+S(1)/S(2)*e**S(3)*(S(4)*c*d-b*e)*x**S(2)/c**S(2)+S(1)/S(3)*e**S(4)*x**S(3)/c+S(1)/S(2)*e*(S(2)*c*d-b*e)*(S(2)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(b*d+a*e))*log(a+b*x+c*x**S(2))/c**S(4)-(S(2)*c**S(4)*d**S(4)+b**S(4)*e**S(4)-S(4)*b**S(2)*c*e**S(3)*(b*d+a*e)-S(4)*c**S(3)*d**S(2)*e*(b*d+S(3)*a*e)+S(2)*c**S(2)*e**S(2)*(S(3)*b**S(2)*d**S(2)+S(6)*a*b*d*e+a**S(2)*e**S(2)))*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c**S(4)*sqrt(b**S(2)-S(4)*a*c))],
[(d+e*x)**S(3)/(a+b*x+c*x**S(2)),x,S(7),e**S(2)*(S(3)*c*d-b*e)*x/c**S(2)+S(1)/S(2)*e**S(3)*x**S(2)/c+S(1)/S(2)*e*(S(3)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(3)*b*d+a*e))*log(a+b*x+c*x**S(2))/c**S(3)-(S(2)*c*d-b*e)*(c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(b*d+S(3)*a*e))*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c**S(3)*sqrt(b**S(2)-S(4)*a*c))],
[(d+e*x)**S(2)/(a+b*x+c*x**S(2)),x,S(7),e**S(2)*x/c+S(1)/S(2)*e*(S(2)*c*d-b*e)*log(a+b*x+c*x**S(2))/c**S(2)-(S(2)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(b*d+a*e))*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c**S(2)*sqrt(b**S(2)-S(4)*a*c))],
[(d+e*x)/(a+b*x+c*x**S(2)),x,S(5),S(1)/S(2)*e*log(a+b*x+c*x**S(2))/c-(S(2)*c*d-b*e)*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/(c*sqrt(b**S(2)-S(4)*a*c))],
[S(1)/(a+b*x+c*x**S(2)),x,S(2),-S(2)*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c)],
[S(1)/((d+e*x)*(a+b*x+c*x**S(2))),x,S(7),e*log(d+e*x)/(c*d**S(2)-b*d*e+a*e**S(2))-S(1)/S(2)*e*log(a+b*x+c*x**S(2))/(c*d**S(2)-b*d*e+a*e**S(2))-(S(2)*c*d-b*e)*arctanh((b+S(2)*c*x)/sqrt(b**S(2)-S(4)*a*c))/((c*d**S(2)-b*d*e+a*e**S(2))*sqrt(b**S(2)-S(4)*a*c))],

# Integrands of the form (d+e x)**(m/S(2)) (a+b x+c x**S(2))**p

# p>S(0)
[(d+e*x)**(S(5)/S(2))*(a+b*x+c*x**S(2)),x,S(2),S(2)/S(7)*(c*d**S(2)-b*d*e+a*e**S(2))*(d+e*x)**(S(7)/S(2))/e**S(3)-S(2)/S(9)*(S(2)*c*d-b*e)*(d+e*x)**(S(9)/S(2))/e**S(3)+S(2)/S(11)*c*(d+e*x)**(S(11)/S(2))/e**S(3)],
[(d+e*x)**(S(3)/S(2))*(a+b*x+c*x**S(2)),x,S(2),S(2)/S(5)*(c*d**S(2)-b*d*e+a*e**S(2))*(d+e*x)**(S(5)/S(2))/e**S(3)-S(2)/S(7)*(S(2)*c*d-b*e)*(d+e*x)**(S(7)/S(2))/e**S(3)+S(2)/S(9)*c*(d+e*x)**(S(9)/S(2))/e**S(3)],

# p<S(0)
[(d+e*x)**(S(5)/S(2))/(a+b*x+c*x**S(2)),x,S(6),S(2)/S(3)*e*(d+e*x)**(S(3)/S(2))/c+S(2)*e*(S(2)*c*d-b*e)*sqrt(d+e*x)/c**S(2)-arctanh(sqrt(S(2))*sqrt(c)*sqrt(d+e*x)/sqrt(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))))*sqrt(S(2))*(S(2)*c**S(3)*d**S(3)-b**S(2)*e**S(3)*(b-sqrt(b**S(2)-S(4)*a*c))-S(3)*c**S(2)*d*e*(b*d+S(2)*a*e-d*sqrt(b**S(2)-S(4)*a*c))+c*e**S(2)*(S(3)*b**S(2)*d+S(3)*a*b*e-S(3)*b*d*sqrt(b**S(2)-S(4)*a*c)-a*e*sqrt(b**S(2)-S(4)*a*c)))/(c**(S(5)/S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))))+arctanh(sqrt(S(2))*sqrt(c)*sqrt(d+e*x)/sqrt(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))*sqrt(S(2))*(S(2)*c**S(3)*d**S(3)-b**S(2)*e**S(3)*(b+sqrt(b**S(2)-S(4)*a*c))-S(3)*c**S(2)*d*e*(b*d+S(2)*a*e+d*sqrt(b**S(2)-S(4)*a*c))+c*e**S(2)*(S(3)*b**S(2)*d+a*e*sqrt(b**S(2)-S(4)*a*c)+S(3)*b*(a*e+d*sqrt(b**S(2)-S(4)*a*c))))/(c**(S(5)/S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))],
[(d+e*x)**(S(3)/S(2))/(a+b*x+c*x**S(2)),x,S(5),S(2)*e*sqrt(d+e*x)/c-arctanh(sqrt(S(2))*sqrt(c)*sqrt(d+e*x)/sqrt(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))))*sqrt(S(2))*(S(2)*c**S(2)*d**S(2)+b*e**S(2)*(b-sqrt(b**S(2)-S(4)*a*c))-S(2)*c*e*(b*d+a*e-d*sqrt(b**S(2)-S(4)*a*c)))/(c**(S(3)/S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))))+arctanh(sqrt(S(2))*sqrt(c)*sqrt(d+e*x)/sqrt(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))*sqrt(S(2))*(S(2)*c**S(2)*d**S(2)+b*e**S(2)*(b+sqrt(b**S(2)-S(4)*a*c))-S(2)*c*e*(b*d+a*e+d*sqrt(b**S(2)-S(4)*a*c)))/(c**(S(3)/S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))],

#  b**S(2)-S(4)*a*c<S(0) 
[(S(1)+S(2)*x)**(S(7)/S(2))/(S(2)+S(3)*x+S(5)*x**S(2)),x,S(15),S(16)/S(75)*(S(1)+S(2)*x)**(S(3)/S(2))+S(4)/S(25)*(S(1)+S(2)*x)**(S(5)/S(2))-S(76)/S(125)*sqrt(S(1)+S(2)*x)+S(1)/S(125)*arctan(sqrt(S(5)/S(2)/(-S(2)+sqrt(S(35))))*(-S(2)*sqrt(S(1)+S(2)*x)+sqrt(S(2)/S(5)*(S(2)+sqrt(S(35))))))*sqrt(S(2)/S(155)*(-S(168698)+S(42875)*sqrt(S(35))))-S(1)/S(125)*arctan(sqrt(S(5)/S(2)/(-S(2)+sqrt(S(35))))*(S(2)*sqrt(S(1)+S(2)*x)+sqrt(S(2)/S(5)*(S(2)+sqrt(S(35))))))*sqrt(S(2)/S(155)*(-S(168698)+S(42875)*sqrt(S(35))))-S(1)/S(125)*log(S(5)*(S(1)+S(2)*x)+sqrt(S(35))-sqrt(S(1)+S(2)*x)*sqrt(S(10)*(S(2)+sqrt(S(35)))))*sqrt(S(1)/S(310)*(S(168698)+S(42875)*sqrt(S(35))))+S(1)/S(125)*log(S(5)*(S(1)+S(2)*x)+sqrt(S(35))+sqrt(S(1)+S(2)*x)*sqrt(S(10)*(S(2)+sqrt(S(35)))))*sqrt(S(1)/S(310)*(S(168698)+S(42875)*sqrt(S(35))))],
[(S(1)+S(2)*x)**(S(5)/S(2))/(S(2)+S(3)*x+S(5)*x**S(2)),x,S(14),S(4)/S(15)*(S(1)+S(2)*x)**(S(3)/S(2))+S(16)/S(25)*sqrt(S(1)+S(2)*x)+S(1)/S(25)*log(S(5)*(S(1)+S(2)*x)+sqrt(S(35))-sqrt(S(1)+S(2)*x)*sqrt(S(10)*(S(2)+sqrt(S(35)))))*sqrt(S(1)/S(310)*(-S(7162)+S(1225)*sqrt(S(35))))-S(1)/S(25)*log(S(5)*(S(1)+S(2)*x)+sqrt(S(35))+sqrt(S(1)+S(2)*x)*sqrt(S(10)*(S(2)+sqrt(S(35)))))*sqrt(S(1)/S(310)*(-S(7162)+S(1225)*sqrt(S(35))))+S(1)/S(25)*arctan(sqrt(S(5)/S(2)/(-S(2)+sqrt(S(35))))*(-S(2)*sqrt(S(1)+S(2)*x)+sqrt(S(2)/S(5)*(S(2)+sqrt(S(35))))))*sqrt(S(2)/S(155)*(S(7162)+S(1225)*sqrt(S(35))))-S(1)/S(25)*arctan(sqrt(S(5)/S(2)/(-S(2)+sqrt(S(35))))*(S(2)*sqrt(S(1)+S(2)*x)+sqrt(S(2)/S(5)*(S(2)+sqrt(S(35))))))*sqrt(S(2)/S(155)*(S(7162)+S(1225)*sqrt(S(35))))],

# Integrands of the form (d+e x)**(m/S(3)) (a+b x+c x**S(2))**p
[(S(3)-x+x**S(2))/x**(S(1)/S(3)),x,S(2),S(9)/S(2)*x**(S(2)/S(3))-S(3)/S(5)*x**(S(5)/S(3))+S(3)/S(8)*x**(S(8)/S(3))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**S(3)*sqrt(a+b*x+c*x**S(2)),x,S(5),S(1)/S(5)*e*(d+e*x)**S(2)*(a+b*x+c*x**S(2))**(S(3)/S(2))/c+S(1)/S(240)*e*(S(192)*c**S(2)*d**S(2)+S(35)*b**S(2)*e**S(2)-S(2)*c*e*(S(75)*b*d+S(16)*a*e)+S(42)*c*e*(S(2)*c*d-b*e)*x)*(a+b*x+c*x**S(2))**(S(3)/S(2))/c**S(3)-S(1)/S(256)*(b**S(2)-S(4)*a*c)*(S(2)*c*d-b*e)*(S(16)*c**S(2)*d**S(2)+S(7)*b**S(2)*e**S(2)-S(4)*c*e*(S(4)*b*d+S(3)*a*e))*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/c**(S(9)/S(2))+S(1)/S(128)*(S(2)*c*d-b*e)*(S(16)*c**S(2)*d**S(2)+S(7)*b**S(2)*e**S(2)-S(4)*c*e*(S(4)*b*d+S(3)*a*e))*(b+S(2)*c*x)*sqrt(a+b*x+c*x**S(2))/c**S(4)],
[(d+e*x)**S(2)*sqrt(a+b*x+c*x**S(2)),x,S(5),S(5)/S(24)*e*(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(3)/S(2))/c**S(2)+S(1)/S(4)*e*(d+e*x)*(a+b*x+c*x**S(2))**(S(3)/S(2))/c-S(1)/S(128)*(b**S(2)-S(4)*a*c)*(S(16)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-S(4)*c*e*(S(4)*b*d+a*e))*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/c**(S(7)/S(2))+S(1)/S(64)*(S(16)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-S(4)*c*e*(S(4)*b*d+a*e))*(b+S(2)*c*x)*sqrt(a+b*x+c*x**S(2))/c**S(3)],

# p<S(0)
[(d+e*x)**S(4)/sqrt(a+b*x+c*x**S(2)),x,S(5),S(1)/S(128)*(S(128)*c**S(4)*d**S(4)+S(35)*b**S(4)*e**S(4)-S(128)*c**S(3)*d**S(2)*e*(S(2)*b*d+S(3)*a*e)-S(40)*b**S(2)*c*e**S(3)*(S(4)*b*d+S(3)*a*e)+S(48)*c**S(2)*e**S(2)*(S(6)*b**S(2)*d**S(2)+S(8)*a*b*d*e+a**S(2)*e**S(2)))*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/c**(S(9)/S(2))+S(7)/S(24)*e*(S(2)*c*d-b*e)*(d+e*x)**S(2)*sqrt(a+b*x+c*x**S(2))/c**S(2)+S(1)/S(4)*e*(d+e*x)**S(3)*sqrt(a+b*x+c*x**S(2))/c+S(1)/S(192)*e*(S(608)*c**S(3)*d**S(3)-S(105)*b**S(3)*e**S(3)+S(20)*b*c*e**S(2)*(S(24)*b*d+S(11)*a*e)-S(8)*c**S(2)*d*e*(S(101)*b*d+S(64)*a*e)+S(2)*c*e*(S(104)*c**S(2)*d**S(2)+S(35)*b**S(2)*e**S(2)-S(4)*c*e*(S(26)*b*d+S(9)*a*e))*x)*sqrt(a+b*x+c*x**S(2))/c**S(4)],
[(d+e*x)**S(3)/sqrt(a+b*x+c*x**S(2)),x,S(4),S(1)/S(16)*(S(2)*c*d-b*e)*(S(8)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-S(4)*c*e*(S(2)*b*d+S(3)*a*e))*arctanh(S(1)/S(2)*(b+S(2)*c*x)/(sqrt(c)*sqrt(a+b*x+c*x**S(2))))/c**(S(7)/S(2))+S(1)/S(3)*e*(d+e*x)**S(2)*sqrt(a+b*x+c*x**S(2))/c+S(1)/S(24)*e*(S(64)*c**S(2)*d**S(2)+S(15)*b**S(2)*e**S(2)-S(2)*c*e*(S(27)*b*d+S(8)*a*e)+S(10)*c*e*(S(2)*c*d-b*e)*x)*sqrt(a+b*x+c*x**S(2))/c**S(3)],

# Integrands of the form (d+e x)**(m/S(2)) (a+b x+c x**S(2))**(p/S(2))

# p>S(0)
[(d+e*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2)),x,S(7),S(2)/S(7)*e*(a+b*x+c*x**S(2))**(S(3)/S(2))*sqrt(d+e*x)/c+S(2)/S(105)*(S(3)*c**S(2)*d**S(2)-S(4)*b**S(2)*e**S(2)+c*e*(S(9)*b*d-S(5)*a*e)+S(12)*c*e*(S(2)*c*d-b*e)*x)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2))/(c**S(2)*e)-S(1)/S(105)*(S(2)*c*d-b*e)*(S(3)*c**S(2)*d**S(2)+S(8)*b**S(2)*e**S(2)-c*e*(S(3)*b*d+S(29)*a*e))*EllipticE(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(d+e*x)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(3)*e**S(2)*sqrt(a+b*x+c*x**S(2))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))+S(4)/S(105)*(c*d**S(2)-b*d*e+a*e**S(2))*(S(3)*c**S(2)*d**S(2)+S(2)*b**S(2)*e**S(2)-c*e*(S(3)*b*d+S(5)*a*e))*EllipticF(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))/(c**S(3)*e**S(2)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2)))],
[(d+e*x)**(S(1)/S(2))*sqrt(a+b*x+c*x**S(2)),x,S(7),S(2)/S(5)*(d+e*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2))/e-S(2)/S(15)*(S(2)*c*d-b*e)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2))/(c*e)-S(2)/S(15)*(c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(b*d+S(3)*a*e))*EllipticE(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(d+e*x)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(2)*e**S(2)*sqrt(a+b*x+c*x**S(2))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))+S(2)/S(15)*(S(2)*c*d-b*e)*(c*d**S(2)-b*d*e+a*e**S(2))*EllipticF(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))/(c**S(2)*e**S(2)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2)))],

# p<S(0)
[(d+e*x)**(S(7)/S(2))/sqrt(a+b*x+c*x**S(2)),x,S(8),S(12)/S(35)*e*(S(2)*c*d-b*e)*(d+e*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2))/c**S(2)+S(2)/S(7)*e*(d+e*x)**(S(5)/S(2))*sqrt(a+b*x+c*x**S(2))/c+S(2)/S(105)*e*(S(71)*c**S(2)*d**S(2)+S(24)*b**S(2)*e**S(2)-c*e*(S(71)*b*d+S(25)*a*e))*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2))/c**S(3)+S(8)/S(105)*(S(2)*c*d-b*e)*(S(11)*c**S(2)*d**S(2)+S(6)*b**S(2)*e**S(2)-c*e*(S(11)*b*d+S(13)*a*e))*EllipticE(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(d+e*x)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(4)*sqrt(a+b*x+c*x**S(2))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))-S(2)/S(105)*(c*d**S(2)-b*d*e+a*e**S(2))*(S(71)*c**S(2)*d**S(2)+S(24)*b**S(2)*e**S(2)-c*e*(S(71)*b*d+S(25)*a*e))*EllipticF(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))/(c**S(4)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2)))],
[(d+e*x)**(S(5)/S(2))/sqrt(a+b*x+c*x**S(2)),x,S(7),S(2)/S(5)*e*(d+e*x)**(S(3)/S(2))*sqrt(a+b*x+c*x**S(2))/c+S(8)/S(15)*e*(S(2)*c*d-b*e)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2))/c**S(2)+S(1)/S(15)*(S(23)*c**S(2)*d**S(2)+S(8)*b**S(2)*e**S(2)-c*e*(S(23)*b*d+S(9)*a*e))*EllipticE(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(d+e*x)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))/(c**S(3)*sqrt(a+b*x+c*x**S(2))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))-S(8)/S(15)*(S(2)*c*d-b*e)*(c*d**S(2)-b*d*e+a*e**S(2))*EllipticF(sqrt((b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/sqrt(b**S(2)-S(4)*a*c))/sqrt(S(2)),sqrt(-S(2)*e*sqrt(b**S(2)-S(4)*a*c)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c)*sqrt(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))*sqrt(c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))/(c**S(3)*sqrt(d+e*x)*sqrt(a+b*x+c*x**S(2)))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**(p/S(3))

# p>S(0)
[(d+e*x)**S(2)*(a+b*x+c*x**S(2))**(S(4)/S(3)),x,S(6),S(15)/S(119)*e*(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(7)/S(3))/c**S(2)+S(3)/S(17)*e*(d+e*x)*(a+b*x+c*x**S(2))**(S(7)/S(3))/c-S(3)/S(935)*(b**S(2)-S(4)*a*c)*(S(17)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-c*e*(S(17)*b*d+S(3)*a*e))*(a+b*x+c*x**S(2))**(S(1)/S(3))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**S(4)*(b+S(2)*c*x))+S(3)/S(374)*(S(17)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-c*e*(S(17)*b*d+S(3)*a*e))*(a+b*x+c*x**S(2))**(S(4)/S(3))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**S(3)*(b+S(2)*c*x))+S(1)/S(935)*S(2)**(S(1)/S(3))*S(3)**(S(3)/S(4))*(b**S(2)-S(4)*a*c)**S(2)*(S(17)*c**S(2)*d**S(2)+S(5)*b**S(2)*e**S(2)-c*e*(S(17)*b*d+S(3)*a*e))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticF((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(S(2)+sqrt(S(3)))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(c**(S(13)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))],
[(d+e*x)*(a+b*x+c*x**S(2))**(S(4)/S(3)),x,S(5),S(3)/S(14)*e*(a+b*x+c*x**S(2))**(S(7)/S(3))/c-S(3)/S(110)*(b**S(2)-S(4)*a*c)*(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(1)/S(3))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**S(3)*(b+S(2)*c*x))+S(3)/S(44)*(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(4)/S(3))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**S(2)*(b+S(2)*c*x))+S(1)/S(55)*S(3)**(S(3)/S(4))*(b**S(2)-S(4)*a*c)**S(2)*(S(2)*c*d-b*e)*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticF((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(S(2)+sqrt(S(3)))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(S(2)**(S(2)/S(3))*c**(S(10)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))],

# p<S(0)
[(d+e*x)**S(3)/(a+b*x+c*x**S(2))**(S(7)/S(3)),x,S(6),-S(3)/S(4)*(d+e*x)**S(2)*(b*d-S(2)*a*e+(S(2)*c*d-b*e)*x)/((b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**(S(4)/S(3)))+S(3)/S(4)*(S(10)*b*c*d*(c*d**S(2)+S(3)*a*e**S(2))-S(8)*a*c*e*(S(2)*c*d**S(2)+S(3)*a*e**S(2))-b**S(2)*(S(11)*c*d**S(2)*e-a*e**S(3))+(S(2)*c*d-b*e)*(S(10)*c**S(2)*d**S(2)-b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(7)*a*e))*x)/(c*(b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2))**(S(1)/S(3)))-S(3)/S(2)*(S(2)*c*d-b*e)*(S(5)*c**S(2)*d**S(2)-b**S(2)*e**S(2)-c*e*(S(5)*b*d-S(9)*a*e))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(S(2)**(S(1)/S(3))*c**(S(5)/S(3))*(b**S(2)-S(4)*a*c)**S(2)*(b+S(2)*c*x)*(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))))-S(3)**(S(3)/S(4))*(S(2)*c*d-b*e)*(S(5)*c**S(2)*d**S(2)-b**S(2)*e**S(2)-c*e*(S(5)*b*d-S(9)*a*e))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticF((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(S(2)**(S(5)/S(6))*c**(S(5)/S(3))*(b**S(2)-S(4)*a*c)**(S(5)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))+S(3)/S(4)*S(3)**(S(1)/S(4))*(S(2)*c*d-b*e)*(S(5)*c**S(2)*d**S(2)-b**S(2)*e**S(2)-c*e*(S(5)*b*d-S(9)*a*e))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticE((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(S(2)-sqrt(S(3)))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(S(2)**(S(1)/S(3))*c**(S(5)/S(3))*(b**S(2)-S(4)*a*c)**(S(5)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))],
[(d+e*x)**S(2)/(a+b*x+c*x**S(2))**(S(7)/S(3)),x,S(6),-S(3)/S(4)*(d+e*x)*(b*d-S(2)*a*e+(S(2)*c*d-b*e)*x)/((b**S(2)-S(4)*a*c)*(a+b*x+c*x**S(2))**(S(4)/S(3)))-S(3)/S(2)*(S(4)*b**S(2)*d*e+S(4)*a*c*d*e-S(5)*b*(c*d**S(2)+a*e**S(2))-(S(10)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(3)*a*e))*x)/((b**S(2)-S(4)*a*c)**S(2)*(a+b*x+c*x**S(2))**(S(1)/S(3)))-S(3)/S(2)*(S(10)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(3)*a*e))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(b**S(2)-S(4)*a*c)**S(2)*(b+S(2)*c*x)*(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))))-S(3)**(S(3)/S(4))*(S(10)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(3)*a*e))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticF((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(S(2)**(S(5)/S(6))*c**(S(2)/S(3))*(b**S(2)-S(4)*a*c)**(S(5)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))+S(3)/S(4)*S(3)**(S(1)/S(4))*(S(10)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(3)*a*e))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))*EllipticE((S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)-sqrt(S(3))))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3)))),sqrt(-S(7)-S(4)*sqrt(S(3))))*sqrt((b+S(2)*c*x)**S(2))*sqrt(S(2)-sqrt(S(3)))*sqrt(((b**S(2)-S(4)*a*c)**(S(2)/S(3))-S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+S(2)*S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(a+b*x+c*x**S(2))**(S(2)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2))/(S(2)**(S(1)/S(3))*c**(S(2)/S(3))*(b**S(2)-S(4)*a*c)**(S(5)/S(3))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))*sqrt((b**S(2)-S(4)*a*c)**(S(1)/S(3))*((b**S(2)-S(4)*a*c)**(S(1)/S(3))+S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3)))/(S(2)**(S(2)/S(3))*c**(S(1)/S(3))*(a+b*x+c*x**S(2))**(S(1)/S(3))+(b**S(2)-S(4)*a*c)**(S(1)/S(3))*(S(1)+sqrt(S(3))))**S(2)))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**(p/S(4))

# p>S(0)
[(d+e*x)*(a+b*x+c*x**S(2))**(S(1)/S(4)),x,S(4),S(1)/S(6)*(S(2)*c*d-b*e)*(b+S(2)*c*x)*(a+b*x+c*x**S(2))**(S(1)/S(4))/c**S(2)+S(2)/S(5)*e*(a+b*x+c*x**S(2))**(S(5)/S(4))/c-S(1)/S(12)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*(S(2)*c*d-b*e)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(9)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],
[(a+b*x+c*x**S(2))**(S(1)/S(4)),x,S(3),S(1)/S(3)*(b+S(2)*c*x)*(a+b*x+c*x**S(2))**(S(1)/S(4))/c-S(1)/S(6)*(b**S(2)-S(4)*a*c)**(S(5)/S(4))*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(5)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],
[(a+b*x+c*x**S(2))**(S(1)/S(4))/(d+e*x),x,S(19),S(2)*(a+b*x+c*x**S(2))**(S(1)/S(4))/e-(-b**S(2)+S(4)*a*c)**(S(3)/S(4))*(c*d**S(2)-b*d*e+a*e**S(2))**(S(1)/S(4))*(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))**(S(3)/S(4))*arctan((-b**S(2)+S(4)*a*c)**(S(1)/S(4))*(S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**(S(1)/S(4))*sqrt(e)/(c**(S(1)/S(4))*(c*d**S(2)-b*d*e+a*e**S(2))**(S(1)/S(4))*sqrt(S(2))))/(c**(S(3)/S(4))*e**(S(3)/S(2))*(a+b*x+c*x**S(2))**(S(3)/S(4)))-(-b**S(2)+S(4)*a*c)**(S(3)/S(4))*(c*d**S(2)-b*d*e+a*e**S(2))**(S(1)/S(4))*(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))**(S(3)/S(4))*arctanh((-b**S(2)+S(4)*a*c)**(S(1)/S(4))*(S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**(S(1)/S(4))*sqrt(e)/(c**(S(1)/S(4))*(c*d**S(2)-b*d*e+a*e**S(2))**(S(1)/S(4))*sqrt(S(2))))/(c**(S(3)/S(4))*e**(S(3)/S(2))*(a+b*x+c*x**S(2))**(S(3)/S(4)))-(b**S(2)-S(4)*a*c)*(S(2)*c*d-b*e)*(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))**(S(3)/S(4))*EllipticPi((S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**(S(1)/S(4)),-S(1)/S(2)*e*sqrt(-b**S(2)+S(4)*a*c)/(sqrt(c)*sqrt(c*d**S(2)-b*d*e+a*e**S(2))),I)*sqrt((b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/(c*e**S(2)*(b+S(2)*c*x)*(a+b*x+c*x**S(2))**(S(3)/S(4))*sqrt(S(2)))-(b**S(2)-S(4)*a*c)*(S(2)*c*d-b*e)*(-c*(a+b*x+c*x**S(2))/(b**S(2)-S(4)*a*c))**(S(3)/S(4))*EllipticPi((S(1)-(b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))**(S(1)/S(4)),S(1)/S(2)*e*sqrt(-b**S(2)+S(4)*a*c)/(sqrt(c)*sqrt(c*d**S(2)-b*d*e+a*e**S(2))),I)*sqrt((b+S(2)*c*x)**S(2)/(b**S(2)-S(4)*a*c))/(c*e**S(2)*(b+S(2)*c*x)*(a+b*x+c*x**S(2))**(S(3)/S(4))*sqrt(S(2)))-(b**S(2)-S(4)*a*c)**(S(1)/S(4))*(S(2)*c*d-b*e)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(1)/S(4))*e**S(2)*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],

# p<S(0)
[(d+e*x)**S(2)/(a+b*x+c*x**S(2))**(S(1)/S(4)),x,S(6),S(7)/S(15)*e*(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(3)/S(4))/c**S(2)+S(2)/S(5)*e*(d+e*x)*(a+b*x+c*x**S(2))**(S(3)/S(4))/c+S(1)/S(10)*(S(20)*c**S(2)*d**S(2)+S(7)*b**S(2)*e**S(2)-S(4)*c*e*(S(5)*b*d+S(2)*a*e))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**(S(5)/S(2))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c)))-S(1)/S(10)*(b**S(2)-S(4)*a*c)**(S(3)/S(4))*(S(20)*c**S(2)*d**S(2)+S(7)*b**S(2)*e**S(2)-S(4)*c*e*(S(5)*b*d+S(2)*a*e))*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(11)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))+S(1)/S(20)*(b**S(2)-S(4)*a*c)**(S(3)/S(4))*(S(20)*c**S(2)*d**S(2)+S(7)*b**S(2)*e**S(2)-S(4)*c*e*(S(5)*b*d+S(2)*a*e))*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(11)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],
[(d+e*x)/(a+b*x+c*x**S(2))**(S(1)/S(4)),x,S(5),S(2)/S(3)*e*(a+b*x+c*x**S(2))**(S(3)/S(4))/c+(S(2)*c*d-b*e)*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/(c**(S(3)/S(2))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c)))-(b**S(2)-S(4)*a*c)**(S(3)/S(4))*(S(2)*c*d-b*e)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(7)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))+S(1)/S(2)*(b**S(2)-S(4)*a*c)**(S(3)/S(4))*(S(2)*c*d-b*e)*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(7)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],
[S(1)/(a+b*x+c*x**S(2))**(S(1)/S(4)),x,S(4),S(2)*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt((b+S(2)*c*x)**S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b+S(2)*c*x)*sqrt(c)*sqrt(b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c)))+(b**S(2)-S(4)*a*c)**(S(3)/S(4))*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticF(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(3)/S(4))*(b+S(2)*c*x)*sqrt(S(2))*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))-(b**S(2)-S(4)*a*c)**(S(3)/S(4))*sqrt(cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))**S(2))/cos(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4))))*EllipticE(sin(S(2)*arctan(c**(S(1)/S(4))*(a+b*x+c*x**S(2))**(S(1)/S(4))*sqrt(S(2))/(b**S(2)-S(4)*a*c)**(S(1)/S(4)))),sqrt(S(1)/S(2)))*sqrt(S(2))*sqrt((b+S(2)*c*x)**S(2))*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))*sqrt((b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2)))/((b**S(2)-S(4)*a*c)*(S(1)+S(2)*sqrt(c)*sqrt(a+b*x+c*x**S(2))/sqrt(b**S(2)-S(4)*a*c))**S(2)))/(c**(S(3)/S(4))*(b+S(2)*c*x)*sqrt(b**S(2)-S(4)*a*c+S(4)*c*(a+b*x+c*x**S(2))))],

# Integrands of the form (d+e x)**(m/S(2)) (a+b x+c x**S(2))**(p/S(4))
[S(1)/((d+e*x)**(S(3)/S(2))*(a+b*x+c*x**S(2))**(S(1)/S(4))),x,S(1),S(2)*hypergeom([-S(1)/S(2),S(1)/S(4)],[S(1)/S(2)],-S(4)*c*(d+e*x)*sqrt(b**S(2)-S(4)*a*c)/((b+S(2)*c*x-sqrt(b**S(2)-S(4)*a*c))*(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))*(b+S(2)*c*x-sqrt(b**S(2)-S(4)*a*c))*((S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c)))*(b+S(2)*c*x+sqrt(b**S(2)-S(4)*a*c))/((b+S(2)*c*x-sqrt(b**S(2)-S(4)*a*c))*(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c)))))**(S(1)/S(4))/((a+b*x+c*x**S(2))**(S(1)/S(4))*(S(2)*c*d-b*e+e*sqrt(b**S(2)-S(4)*a*c))*sqrt(d+e*x))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p with m symbolic
[(d+e*x)**m*(a+b*x+c*x**S(2))**S(4),x,S(2),(c*d**S(2)-b*d*e+a*e**S(2))**S(4)*(d+e*x)**(S(1)+m)/(e**S(9)*(S(1)+m))-S(4)*(S(2)*c*d-b*e)*(c*d**S(2)-b*d*e+a*e**S(2))**S(3)*(d+e*x)**(S(2)+m)/(e**S(9)*(S(2)+m))+S(2)*(c*d**S(2)-b*d*e+a*e**S(2))**S(2)*(S(14)*c**S(2)*d**S(2)+S(3)*b**S(2)*e**S(2)-S(2)*c*e*(S(7)*b*d-a*e))*(d+e*x)**(S(3)+m)/(e**S(9)*(S(3)+m))-S(4)*(S(2)*c*d-b*e)*(c*d**S(2)-b*d*e+a*e**S(2))*(S(7)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(7)*b*d-S(3)*a*e))*(d+e*x)**(S(4)+m)/(e**S(9)*(S(4)+m))+(S(70)*c**S(4)*d**S(4)+b**S(4)*e**S(4)-S(4)*b**S(2)*c*e**S(3)*(S(5)*b*d-S(3)*a*e)-S(20)*c**S(3)*d**S(2)*e*(S(7)*b*d-S(3)*a*e)+S(6)*c**S(2)*e**S(2)*(S(15)*b**S(2)*d**S(2)-S(10)*a*b*d*e+a**S(2)*e**S(2)))*(d+e*x)**(S(5)+m)/(e**S(9)*(S(5)+m))-S(4)*c*(S(2)*c*d-b*e)*(S(7)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(7)*b*d-S(3)*a*e))*(d+e*x)**(S(6)+m)/(e**S(9)*(S(6)+m))+S(2)*c**S(2)*(S(14)*c**S(2)*d**S(2)+S(3)*b**S(2)*e**S(2)-S(2)*c*e*(S(7)*b*d-a*e))*(d+e*x)**(S(7)+m)/(e**S(9)*(S(7)+m))-S(4)*c**S(3)*(S(2)*c*d-b*e)*(d+e*x)**(S(8)+m)/(e**S(9)*(S(8)+m))+c**S(4)*(d+e*x)**(S(9)+m)/(e**S(9)*(S(9)+m))],
[(d+e*x)**m*(a+b*x+c*x**S(2))**S(3),x,S(2),(c*d**S(2)-b*d*e+a*e**S(2))**S(3)*(d+e*x)**(S(1)+m)/(e**S(7)*(S(1)+m))-S(3)*(S(2)*c*d-b*e)*(c*d**S(2)-b*d*e+a*e**S(2))**S(2)*(d+e*x)**(S(2)+m)/(e**S(7)*(S(2)+m))+S(3)*(c*d**S(2)-b*d*e+a*e**S(2))*(S(5)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(5)*b*d-a*e))*(d+e*x)**(S(3)+m)/(e**S(7)*(S(3)+m))-(S(2)*c*d-b*e)*(S(10)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-S(2)*c*e*(S(5)*b*d-S(3)*a*e))*(d+e*x)**(S(4)+m)/(e**S(7)*(S(4)+m))+S(3)*c*(S(5)*c**S(2)*d**S(2)+b**S(2)*e**S(2)-c*e*(S(5)*b*d-a*e))*(d+e*x)**(S(5)+m)/(e**S(7)*(S(5)+m))-S(3)*c**S(2)*(S(2)*c*d-b*e)*(d+e*x)**(S(6)+m)/(e**S(7)*(S(6)+m))+c**S(3)*(d+e*x)**(S(7)+m)/(e**S(7)*(S(7)+m))],
[(d+e*x)**m/(a+b*x+c*x**S(2))**(S(5)/S(2)),x,S(2),(d+e*x)**(S(1)+m)*AppellFS(1)(S(1)+m,S(5)/S(2),S(5)/S(2),S(2)+m,S(2)*c*(d+e*x)/(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))),S(2)*c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))*(S(1)-S(2)*c*(d+e*x)/(S(2)*c*d-e*(b-sqrt(b**S(2)-S(4)*a*c))))**(S(5)/S(2))*(S(1)-S(2)*c*(d+e*x)/(S(2)*c*d-e*(b+sqrt(b**S(2)-S(4)*a*c))))**(S(5)/S(2))/(e*(S(1)+m)*(a+b*x+c*x**S(2))**(S(5)/S(2)))],

# Integrands of the form (d+e x)**m (a+b x+c x**S(2))**p with p symbolic
[S(1)/((b*e-c*e*x)**(S(1)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(1)/S(3))),x,S(9),-S(1)/S(6)*(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3))*log(S(1)+c**S(2)*e**(S(2)/S(3))*x**S(2)/(b**S(3)*e-c**S(3)*e*x**S(3))**(S(2)/S(3))-c*e**(S(1)/S(3))*x/(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3)))/(c*e**(S(1)/S(3))*(b*e-c*e*x)**(S(1)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(1)/S(3)))+S(1)/S(3)*(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3))*log(S(1)+c*e**(S(1)/S(3))*x/(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3)))/(c*e**(S(1)/S(3))*(b*e-c*e*x)**(S(1)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(1)/S(3)))-(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3))*arctan((S(1)-S(2)*c*e**(S(1)/S(3))*x/(b**S(3)*e-c**S(3)*e*x**S(3))**(S(1)/S(3)))/sqrt(S(3)))/(c*e**(S(1)/S(3))*(b*e-c*e*x)**(S(1)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(1)/S(3))*sqrt(S(3)))],
[S(1)/((b*e-c*e*x)**(S(2)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(2)/S(3))),x,S(3),x*(S(1)-c**S(3)*x**S(3)/b**S(3))**(S(2)/S(3))*hypergeom([S(1)/S(3),S(2)/S(3)],[S(4)/S(3)],c**S(3)*x**S(3)/b**S(3))/((b*e-c*e*x)**(S(2)/S(3))*(b**S(2)+b*c*x+c**S(2)*x**S(2))**(S(2)/S(3)))],
[(b*e-c*e*x)**p*(b**S(2)+b*c*x+c**S(2)*x**S(2))**p,x,S(3),x*(b*e-c*e*x)**p*(b**S(2)+b*c*x+c**S(2)*x**S(2))**p*hypergeom([S(1)/S(3),-p],[S(4)/S(3)],c**S(3)*x**S(3)/b**S(3))/(S(1)-c**S(3)*x**S(3)/b**S(3))**p]
