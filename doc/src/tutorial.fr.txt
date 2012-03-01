.. _tutorial:

========
Tutoriel
========

.. role:: input(strong)

Introduction
============

SymPy est une bibliothèque Python pour les mathématiques symboliques.  SymPy
prévoit devenir un système complet de calcul formel ("CAS" en anglais : "Computer
Algebra System") tout en gardant le code aussi simple que possible afin qu'il soit
compréhensible et facilement extensible.  SymPy est entièrement écrit en Python et
ne nécessite aucune bibliothèque externe.

Ce tutoriel donne un aperçu et une introduction à SymPy.  Lisez-le pour vous faire
une idée de ce que SymPy peut faire pour vous (et comment).  Pour en savoir plus,
consultez également le :ref:`SymPy User's Guide <guide>`, la
:ref:`SymPy Modules Reference <module-docs>`, ou le
`code source <https://github.com/sympy/sympy/>`_ directement.

Premiers Pas avec SymPy
=======================

La façon la plus simple de télécharger SymPy est de se rendre sur
http://code.google.com/p/sympy/ et de télécharger la dernière archive
des Téléchargements Recommandés ("Featured Downloads" en anglais).

.. image:: figures/featured-downloads.png

Décompressez-la :

.. parsed-literal::

    $ :input:`tar xzf sympy-0.5.12.tar.gz`

et essayez-la dans un interpréteur Python :

.. parsed-literal::

    $ :input:`cd sympy-0.5.12`
    $ :input:`python`
    Python 2.4.4 (#2, Jan  3 2008, 13:36:28)
    [GCC 4.2.3 20071123 (prerelease) (Debian 4.2.2-4)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from sympy import Symbol, cos
    >>> x = Symbol("x")
    >>> (1/cos(x)).series(x, 0, 10)
    1 + x**2/2 + 5*x**4/24 + 61*x**6/720 + 277*x**8/8064 + O(x**10)

Vous pouvez utiliser SymPy comme illustré ci-dessus et c'est en effet
la méthode recommandée si vous voulez l'utiliser dans votre programme.
Vous pouvez aussi l'installer en utilisant ``./setup.py install``
comme tout autre module Python, ou simplement installer un paquet
dans votre distribution Linux préférée, par exemple :

.. topic:: Installer SymPy dans Debian

  .. parsed-literal::

    $ :input:`sudo apt-get install python-sympy`
    Reading package lists... Done
    Building dependency tree
    Reading state information... Done
    The following NEW packages will be installed:
      python-sympy
    0 upgraded, 1 newly installed, 0 to remove and 18 not upgraded.
    Need to get 991kB of archives.
    After this operation, 5976kB of additional disk space will be used.
    Get:1 http://ftp.cz.debian.org unstable/main python-sympy 0.5.12-1 [991kB]
    Fetched 991kB in 2s (361kB/s)
    Selecting previously deselected package python-sympy.
    (Reading database ... 232619 files and directories currently installed.)
    Unpacking python-sympy (from .../python-sympy_0.5.12-1_all.deb) ...
    Setting up python-sympy (0.5.12-1) ...


Pour d'autres moyens d'installer SymPy, consultez l'onglet Téléchargements_
("Downloads" en anglais) sur la page web de SymPy.

.. _Téléchargements: http://code.google.com/p/sympy/wiki/DownloadInstallation?tm=2


La Console isympy
-----------------

Pour tester de nouvelles fonctionnalités, ou si vous cherchez à savoir comment
accomplir certaines choses, vous pouvez utiliser notre emballage (wrapper) spécial
enrobant IPython appelé ``isympy`` (situé dans ``bin/isympy`` si vous vous trouvez
dans le répertoire source).  Cette commande lance un shell Python standard qui
importe les modules SymPy les plus courants, les symboles définis x, y, z et
d'autres choses à l'avance :

.. parsed-literal::

    $ :input:`cd sympy`
    $ :input:`./bin/isympy`
    IPython console for SymPy 0.7.1-git (Python 2.7.1) (ground types: gmpy)

    These commands were executed:
    >>> from __future__ import division
    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)

    Documentation can be found at http://www.sympy.org

    In [1]: :input:`(1/cos(x)).series(x, 0, 10)`
    Out[1]:
         2      4       6        8
        x    5*x    61*x    277*x     / 10\
    1 + ── + ──── + ───── + ────── + O\x  /
        2     24     720     8064

.. note::

    Les commandes que vous entrez sont en gras.  Ainsi ce que nous avons fait en 3
    lignes dans un interpréteur Python habituel peut être fait en une ligne avec
    isympy.


Utiliser SymPy comme une calculatrice
-------------------------------------

SymPy possède trois types numériques natifs :
Float (Flottant), Rational (Rationnel) et Integer (Entier).

La classe Rational représente un nombre rationnel en tant que paire de deux
Integers : le numérateur et le dénominateur.  Donc Rational(1,2) représente 1/2,
Rational(5,2) représente 5/2, et ainsi de suite.

::

    >>> from sympy import Rational
    >>> a = Rational(1,2)

    >>> a
    1/2

    >>> a*2
    1

    >>> Rational(2)**50/Rational(10)**50
    1/88817841970012523233890533447265625


Faites attention lorsque vous travaillez avec les nombres entiers et à virgule
flottante de Python, surtout dans des divisions, puisque vous pourriez créer un
nombre Python, et non pas un nombre SymPy.  Un ratio de deux entiers Python
pourrait créer un flottant -- la "vraie division", standard sous Python 3 et le
comportement par défaut de ``isympy`` qui importe l'opération division du module
__future__ ::

    >>> 1/2 #doctest: +SKIP
    0.5

Cependant, dans les versions antérieures à Python 3, l'opération division résulte
par défaut en une divison euclidienne et entraînera :

    >>> 1/2 #doctest: +SKIP
    0

Dans les deux cas, toutefois, vous n'avez pas affaire à un nombre SymPy parce que
Python a créé son propre nombre.  La plupart du temps vous travaillerez fort
probablement avec des nombres rationels, alors prenez garde à utiliser Rational
pour obtenir le résultat en terme d'objets SymPy.  On pourrait trouver pratique
d'assigner Rational à ``R`` ::

    >>> R = Rational
    >>> R(1, 2)
    1/2
    >>> R(1)/2 # R(1) est un Integer SymPy et Integer/int donne un Rational
    1/2

SymPy offre également quelques constantes spéciales, comme e et pi,
qui sont traitées comme des symboles (1+pi ne sera pas évalué
sous forme numérique, mais restera plutôt sous la forme 1+pi),
et qui ont une précision arbitraire ::

    >>> from sympy import pi, E
    >>> pi**2
    pi**2

    >>> pi.evalf()
    3.14159265358979

    >>> (pi + E).evalf()
    5.85987448204884

Comme vous le voyez, evalf convertit l'expression en un nombre à virgule
flottante.

Le symbole ``oo`` est utilisé pour une classe définissant l'infini mathématique ::

    >>> from sympy import oo
    >>> oo > 99999
    True
    >>> oo + 1
    oo

Symboles
--------

Contrairement à d'autres systèmes de calcul formel,
SymPy vous oblige à déclarer les variables symboliques explicitement ::

    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> y = Symbol('y')

À gauche se trouve la variable Python normale à laquelle est assignée une instance
de la classe Symbol de SymPy.  Certains symboles prédefinis (dont les symboles
portant un nom grec) sont disponibles après importation :

    >>> from sympy.abc import x, theta

Plusieurs symboles peuvent également être déclarés en utilisant une seule commande
avec les fonctions ``symbols`` ou ``var``.  Ces fonctions acceptent une chaîne de
texte contenant une liste ou un interval de symboles à déclarer, la variante
``var`` ajoutant également les symboles à la liste de noms courante (namespace) :

    >>> from sympy import symbols, var
    >>> a, b, c = symbols('a,b,c')
    >>> d, e, f = symbols('d:f')
    >>> var('g:h')
    (g, h)
    >>> var('g:2')
    (g0, g1)

Les instances de la classe Symbol travaillent ensemble et sont les fondements
nécessaires pour construire des expressions ::

    >>> x+y+x-y
    2*x

    >>> (x+y)**2
    (x + y)**2

    >>> ((x+y)**2).expand()
    x**2 + 2*x*y + y**2

Elles peuvent être remplacées par d'autres nombres, symboles ou expressions
grâce à l'utilisation de ``subs(old, new)`` ::

    >>> ((x+y)**2).subs(x, 1)
    (y + 1)**2

    >>> ((x+y)**2).subs(x, y)
    4*y**2

    >>> ((x+y)**2).subs(x, 1-y)
    1

Pour le reste du tutoriel, il est présumé que nous avons exécuté ::

    >>> from sympy import init_printing
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)

Ceci donnera un meilleur aspect aux résultats affichées.  Pour en savoir plus sur
l'affichage, consultez la section :ref:`printing-tutorial` ci-dessous.  Si vous
avez une police unicode installée, vous pouvez passer use_unicode=True pour un
affichage un peu plus agréable.

Algèbre
=======

Pour la décomposition en éléments simples des fractions, utilisez
``apart(expr, x)`` ::

    >>> from sympy import apart
    >>> from sympy.abc import x, y, z

    >>> 1/( (x+2)*(x+1) )
           1
    ---------------
    (x + 1)*(x + 2)

    >>> apart(1/( (x+2)*(x+1) ), x)
        1       1
    - ----- + -----
      x + 2   x + 1

    >>> (x+1)/(x-1)
    x + 1
    -----
    x - 1

    >>> apart((x+1)/(x-1), x)
          2
    1 + -----
        x - 1

Pour remettre tout ensemble, utilisez ``together(expr, x)`` ::

    >>> from sympy import together
    >>> together(1/x + 1/y + 1/z)
    x*y + x*z + y*z
    ---------------
         x*y*z

    >>> together(apart((x+1)/(x-1), x), x)
    x + 1
    -----
    x - 1

    >>> together(apart(1/( (x+2)*(x+1) ), x), x)
           1
    ---------------
    (x + 1)*(x + 2)


.. index:: calculus

Calcul
======

.. index:: limits

Limites
-------

Les limites sont simples à utiliser dans SymPy, elles suivent la syntaxe
``limit(function, variable, point)``, donc pour calculer la limite de f(x) pour x
tendant vers 0, il suffit d'entrer ``limit(f, x, 0)`` ::

   >>> from sympy import limit, Symbol, sin, oo
   >>> x = Symbol("x")
   >>> limit(sin(x)/x, x, 0)
   1

Vous pouvez aussi calculer la limite tendant vers l'infini ::

   >>> limit(x, x, oo)
   oo

   >>> limit(1/x, x, oo)
   0

   >>> limit(x**x, x, 0)
   1

Pour des exemples non-triviaux de limites, vous pouvez lire le fichier de test
`test_demidovich.py
<https://github.com/sympy/sympy/blob/master/sympy/series/tests/test_demidovich.py>`_

.. index:: differentiation, diff

Dérivation
----------

Vous pouvez dériver n'importe quelle expression SymPy en utilisant
``diff(func, var)``.  Exemples ::

    >>> from sympy import diff, Symbol, sin, tan
    >>> x = Symbol('x')
    >>> diff(sin(x), x)
    cos(x)
    >>> diff(sin(2*x), x)
    2*cos(2*x)

    >>> diff(tan(x), x)
       2
    tan (x) + 1

Vous pouvez vérifier que c'est correct avec ::

    >>> from sympy import limit
    >>> from sympy.abc import delta
    >>> limit((tan(x + delta) - tan(x))/delta, delta, 0)
       2
    tan (x) + 1

Des dérivées d'ordre supérieur peuvent être calculées en utilisant la méthode
``diff(func, var, n)`` ::

    >>> diff(sin(2*x), x, 1)
    2*cos(2*x)

    >>> diff(sin(2*x), x, 2)
    -4*sin(2*x)

    >>> diff(sin(2*x), x, 3)
    -8*cos(2*x)


.. index::
    single: series expansion
    single: expansion; series

Développement en série
----------------------

Utilisez ``.series(var, point, order)`` ::

    >>> from sympy import Symbol, cos
    >>> x = Symbol('x')
    >>> cos(x).series(x, 0, 10)
         2    4     6      8
        x    x     x      x      / 10\
    1 - -- + -- - --- + ----- + O\x  /
        2    24   720   40320
    >>> (1/cos(x)).series(x, 0, 10)
         2      4       6        8
        x    5*x    61*x    277*x     / 10\
    1 + -- + ---- + ----- + ------ + O\x  /
        2     24     720     8064

Un autre exemple simple ::

    >>> from sympy import Integral, pprint

    >>> y = Symbol("y")
    >>> e = 1/(x + y)
    >>> s = e.series(x, 0, 5)

    >>> print(s)
    1/y - x/y**2 + x**2/y**3 - x**3/y**4 + x**4/y**5 + O(x**5)
    >>> pprint(s)
              2    3    4
    1   x    x    x    x     / 5\
    - - -- + -- - -- + -- + O\x /
    y    2    3    4    5
        y    y    y    y

.. index:: integration

Intégration
-----------

SymPy supporte l'intégration définie et indéfinie de fonctions
élémentaires et spéciales avec ``integrate()``, qui utilise le
puissant algorithme de Risch-Norman étendu et quelques heuristiques
ainsi que le filtrage par motif ::

    >>> from sympy import integrate, erf, exp, sin, log, oo, pi, sinh, symbols
    >>> x, y = symbols('x,y')

Vous pouvez intégrer des fonctions élémentaires ::

    >>> integrate(6*x**5, x)
     6
    x
    >>> integrate(sin(x), x)
    -cos(x)
    >>> integrate(log(x), x)
    x*log(x) - x
    >>> integrate(2*x + sinh(x), x)
     2
    x  + cosh(x)

Les fonctions spéciales sont aussi facilement gérées ::

    >>> integrate(exp(-x**2)*erf(x), x)
      ____    2
    \/ pi *erf (x)
    --------------
          4

Il est possible de calculer l'intégrale définie ::

    >>> integrate(x**3, (x, -1, 1))
    0
    >>> integrate(sin(x), (x, 0, pi/2))
    1
    >>> integrate(cos(x), (x, -pi/2, pi/2))
    2

Les intégrales impropres sont aussi supportées ::

    >>> integrate(exp(-x), (x, 0, oo))
    1
    >>> integrate(log(x), (x, 0, 1))
    -1

.. index::
    single: complex numbers
    single: expansion; complex

Nombres complexes
-----------------

Outre l'unité imaginaire, I, qui est imaginaire, les symboles peuvent être créés
avec certains attributs (par exemple réel, positif, complexe, etc.) affectant leur
comportement ::

    >>> from sympy import Symbol, exp, I
    >>> x = Symbol("x") # un simple x sans attribut
    >>> exp(I*x).expand()
     I*x
    e
    >>> exp(I*x).expand(complex=True)
       -im(x)               -im(x)
    I*e      *sin(re(x)) + e      *cos(re(x))
    >>> x = Symbol("x", real=True)
    >>> exp(I*x).expand(complex=True)
    I*sin(x) + cos(x)

Fonctions
---------

**trigonométrie** ::

    >>> from sympy import asin, asinh, cos, sin, sinh, symbols, I
    >>> x, y = symbols('x,y')

    >>> sin(x+y).expand(trig=True)
    sin(x)*cos(y) + sin(y)*cos(x)

    >>> cos(x+y).expand(trig=True)
    -sin(x)*sin(y) + cos(x)*cos(y)

    >>> sin(I*x)
    I*sinh(x)

    >>> sinh(I*x)
    I*sin(x)

    >>> asinh(I)
    I*pi
    ----
     2

    >>> asinh(I*x)
    I*asin(x)

    >>> sin(x).series(x, 0, 10)
         3     5     7       9
        x     x     x       x       / 10\
    x - -- + --- - ---- + ------ + O\x  /
        6    120   5040   362880

    >>> sinh(x).series(x, 0, 10)
         3     5     7       9
        x     x     x       x       / 10\
    x + -- + --- + ---- + ------ + O\x  /
        6    120   5040   362880

    >>> asin(x).series(x, 0, 10)
         3      5      7       9
        x    3*x    5*x    35*x     / 10\
    x + -- + ---- + ---- + ----- + O\x  /
        6     40    112     1152

    >>> asinh(x).series(x, 0, 10)
         3      5      7       9
        x    3*x    5*x    35*x     / 10\
    x - -- + ---- - ---- + ----- + O\x  /
        6     40    112     1152

**harmoniques sphériques** ::

    >>> from sympy import Ylm
    >>> from sympy.abc import theta, phi

    >>> Ylm(1, 0, theta, phi)
      ___
    \/ 3 *cos(theta)
    ----------------
            ____
        2*\/ pi

    >>> Ylm(1, 1, theta, phi)
       ___  I*phi
    -\/ 6 *e     *sin(theta)
    ------------------------
                ____
            4*\/ pi

    >>> Ylm(2, 1, theta, phi)
       ____  I*phi
    -\/ 30 *e     *sin(theta)*cos(theta)
    ------------------------------------
                      ____
                  4*\/ pi

**factorielles et fonction gamma** ::

    >>> from sympy import factorial, gamma, Symbol
    >>> x = Symbol("x")
    >>> n = Symbol("n", integer=True)

    >>> factorial(x)
    x!

    >>> factorial(n)
    n!

    >>> gamma(x + 1).series(x, 0, 3) # c'est à dire factorial(x)
                          /          2     2\        
                        2 |EulerGamma    pi |    / 3\
    1 - EulerGamma*x + x *|----------- + ---| + O\x /
                          \     2         12/        

**fonction zeta** ::

    >>> from sympy import zeta
    >>> zeta(4, x)
    zeta(4, x)

    >>> zeta(4, 1)
      4
    pi
    ---
     90

    >>> zeta(4, 2)
           4
         pi
    -1 + ---
          90

    >>> zeta(4, 3)
             4
      17   pi
    - -- + ---
      16    90


**polynômes** ::

    >>> from sympy import assoc_legendre, chebyshevt, legendre, hermite
    >>> chebyshevt(2, x)
       2
    2*x  - 1

    >>> chebyshevt(4, x)
       4      2
    8*x  - 8*x  + 1

    >>> legendre(2, x)
       2
    3*x    1
    ---- - -
     2     2

    >>> legendre(8, x)
          8         6         4        2
    6435*x    3003*x    3465*x    315*x     35
    ------- - ------- + ------- - ------ + ---
      128        32        64       32     128

    >>> assoc_legendre(2, 1, x)
            __________
           /    2
    -3*x*\/  - x  + 1

    >>> assoc_legendre(2, 2, x)
         2
    - 3*x  + 3

    >>> hermite(3, x)
       3
    8*x  - 12*x

.. index:: equations; differential, diff, dsolve

Équations différentielles
-------------------------

Dans ``isympy`` ::

    >>> from sympy import Function, Symbol, dsolve
    >>> f = Function('f')
    >>> x = Symbol('x')
    >>> f(x).diff(x, x) + f(x)
            2
           d
    f(x) + ---(f(x))
             2
           dx

    >>> dsolve(f(x).diff(x, x) + f(x), f(x))
    f(x) = C1*sin(x) + C2*cos(x)

.. index:: equations; algebraic, solve

Équations Algébriques
---------------------

Dans ``isympy`` ::

    >>> from sympy import solve, symbols
    >>> x, y = symbols('x,y')
    >>> solve(x**4 - 1, x)
    [-1, 1, -I, I]

    >>> solve([x + 5*y - 2, -3*x + 6*y - 15], [x, y])
    {x: -3, y: 1}

.. index:: linear algebra

Algèbre Linéaire
================

.. index:: Matrix

Matrices
--------

Les matrices sont créées en instanciant la classe Matrix ::

    >>> from sympy import Matrix, Symbol
    >>> Matrix([[1,0], [0,1]])
    [1  0]
    [    ]
    [0  1]

Elles peuvent aussi contenir des symboles ::

    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> A = Matrix([[1,x], [y,1]])
    >>> A
    [1  x]
    [    ]
    [y  1]

    >>> A**2
    [x*y + 1    2*x  ]
    [                ]
    [  2*y    x*y + 1]

Pour plus d'informations et d'exemples avec les Matrices,
voyez le Tutoriel Algèbre Linéaire.

.. index:: pattern matching, match, Wild, WildFunction

Filtrage par motif
==================

Utilisez la méthode ``.match()``, accompagnée de la classe ``Wild``, pour
effectuer un filtrage par motif sur des expressions.  La méthode renvoit un
dictionnaire avec les substitutions requises, comme suit ::

    >>> from sympy import Symbol, Wild
    >>> x = Symbol('x')
    >>> p = Wild('p')
    >>> (5*x**2).match(p*x**2)
    {p: 5}

    >>> q = Wild('q')
    >>> (x**2).match(p*x**q)
    {p: 1, q: 2}

Si le filtrage est infructueux, ``None`` est renvoyé ::

    >>> print (x+1).match(p**x)
    None

On peut aussi utiliser le paramètre ``exclude`` de la classe ``Wild``
pour s'assurer que certaines choses ne figurent pas dans le résultat ::

    >>> p = Wild('p', exclude=[1,x])
    >>> print (x+1).match(x+p) # 1 est exclu
    None
    >>> print (x+1).match(p+1) # x est exclu
    None
    >>> print (x+1).match(x+2+p) # -1 n'est pas exclu
    {p_: -1}

.. _printing-tutorial:

Affichage
=========

Les expressions peuvent être affichées de plusieures façons.

**Affichage standard**

C'est ce que ``str(expression)`` renvoit et ça ressemble à ça :

    >>> from sympy import Integral
    >>> from sympy.abc import x
    >>> print x**2
    x**2
    >>> print 1/x
    1/x
    >>> print Integral(x**2, x)
    Integral(x**2, x)

**Affichage Amélioré**

Ceci est un affichage amélioré en art ASCII produit par la fonction ``pprint`` :

    >>> from sympy import Integral, pprint
    >>> from sympy.abc import x
    >>> pprint(x**2)
     2
    x
    >>> pprint(1/x)
    1
    -
    x
    >>> pprint(Integral(x**2, x))
      /
     |
     |  2
     | x  dx
     |
    /

Si vous avez une police unicode installée, l'affichage amélioré l'utilisera par
défaut.  Vous pouvez outrepasser ce comportement en utilisant l'option
``use_unicode`` :

    >>> pprint(Integral(x**2, x), use_unicode=True)
    ⌠
    ⎮  2
    ⎮ x  dx
    ⌡


Voyez aussi la page wiki `Pretty Printing
<https://github.com/sympy/sympy/wiki/Pretty-Printing>`_
pour plus d'exemples concernant l'affichage amélioré et l'option unicode.

Astuce : pour activer l'affichage amélioré par défaut
dans l'interpréteur Python, utilisez ::

    $ python
    Python 2.5.2 (r252:60911, Jun 25 2008, 17:58:32)
    [GCC 4.3.1] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from sympy import init_printing, var, Integral
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)
    >>> var("x")
    x
    >>> x**3/3
     3
    x
    --
    3
    >>> Integral(x**2, x) #doctest: +NORMALIZE_WHITESPACE
      /
     |
     |  2
     | x  dx
     |
    /

**Python**

    >>> from sympy.printing.python import python
    >>> from sympy import Integral
    >>> from sympy.abc import x
    >>> print python(x**2)
    x = Symbol('x')
    e = x**2
    >>> print python(1/x)
    x = Symbol('x')
    e = 1/x
    >>> print python(Integral(x**2, x))
    x = Symbol('x')
    e = Integral(x**2, x)


**LaTeX**

    >>> from sympy import Integral, latex
    >>> from sympy.abc import x
    >>> latex(x**2)
    x^{2}
    >>> latex(x**2, mode='inline')
    $x^{2}$
    >>> latex(x**2, mode='equation')
    \begin{equation}x^{2}\end{equation}
    >>> latex(x**2, mode='equation*')
    \begin{equation*}x^{2}\end{equation*}
    >>> latex(1/x)
    \frac{1}{x}
    >>> latex(Integral(x**2, x))
    \int x^{2}\, dx

**MathML**

::

    >>> from sympy.printing.mathml import mathml
    >>> from sympy import Integral, latex
    >>> from sympy.abc import x
    >>> print mathml(x**2)
    <apply><power/><ci>x</ci><cn>2</cn></apply>
    >>> print mathml(1/x)
    <apply><power/><ci>x</ci><cn>-1</cn></apply>

**Pyglet**

    >>> from sympy import Integral, preview
    >>> from sympy.abc import x
    >>> preview(Integral(x**2, x)) #doctest:+SKIP

Si pyglet est installé, une fenêtre pyglet apparaîtra avec l'expression rendue par
LaTeX :

.. image:: pics/pngview1.png

Notes
-----

``isympy`` appelle ``pprint`` automatiquement, c'est donc pourquoi
vous voyez un affichage amélioré par défaut.

Notez qu'il y a aussi un module d'affichage disponible, ``sympy.printing``.
D'autres méthodes d'affichage disponibles dans ce module sont :

* ``pretty(expr)``, ``pretty_print(expr)``, ``pprint(expr)``:
Renvoit ou affiche, respectivement, une représentation améliorée de ``expr``.
C'est la même chose que le second niveau de représentation décrit ci-dessus.

* ``latex(expr)``, ``print_latex(expr)``: Renvoit ou affiche, respectivement,
une représentation `LaTeX <http://www.latex-project.org/>`_ de ``expr``

* ``mathml(expr)``, ``print_mathml(expr)``: Renvoit ou affiche, respectivement,
une représentation `MathML <http://www.w3.org/Math/>`_ de ``expr``.

* ``print_gtk(expr)``: Affiche ``expr`` dans `Gtkmathview
<http://helm.cs.unibo.it/mml-widget/>`_, un widget GTK qui affiche du code
MathML.  Le programme `Gtkmathview <http://helm.cs.unibo.it/mml-widget/>`_ est
requis.

Plus de Documentation
=====================

Il est maintenant temps d'en apprendre plus sur SymPy.  Consultez le
:ref:`SymPy User's Guide <guide>` et la
:ref:`SymPy Modules Reference <module-docs>`.

Veillez également à consulter notre `wiki.sympy.org <http://wiki.sympy.org/>`_
public, qui contient beaucoup d'exemples utiles, des tutoriels, des livres de
recettes auxquels nous et nos utilisateurs ont contribué et que nous vous
encourageons à modifier.
