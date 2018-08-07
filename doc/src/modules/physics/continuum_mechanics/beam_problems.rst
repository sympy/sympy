===========================================================
Solving Beam Bending Problems using Singularity Functions
===========================================================

To make this document easier to read, we are going to enable pretty printing.

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(use_unicode=True, wrap_line=False)

Beam
====

A Beam is a structural element that is capable of withstanding load
primarily by resisting against bending. Beams are characterized by
their second moment of area, their length and their elastic modulus.

In Sympy, we can construct a 2D beam objects with the following properties :

- Length
- Elastic Modulus
- Second Moment of Area
- Variable : A symbol that can be used as a variable along the length. By default,
  this is set to ``Symbol(x)``.
- Boundary Conditions
    - bc_slope : Boundary conditions for slope.
    - bc_deflection : Boundary conditions for deflection.
- Load Distribution

We have following methods under the beam class:

- apply_load
- solve_for_reaction_loads
- shear_force
- bending_moment
- slope
- deflection


Examples
========

Let us solve some beam bending problems using this module :

Example 1
---------

A beam of length 9 meters is having a fixed support at the start.
A distributed constant load of 8 kN/m is applied downward from the starting
point till 5 meters away from the start. A clockwise moment of 50 kN-m is
applied at 5 meters away from the start of the beam. A downward point load
of 12 kN is applied at the end.

::

                                     
  \\\\|    8 KN/m                    
  \\\\|_____________                 
  \\\\|| | | | | | |   ⭯ 50 KN-m     12 KN
  \\\\|V V V V V V V   |               |
  \\\\|________________|_______________V
  \\\\|________________|_______________|
  \\\\|                                :
  \\\\|----------------|---------------|
             5.0 m            4.0 m

.. note::

    Since a user is free to choose its own sign convention we are considering
    the upward forces and clockwise bending moment being positive.


>>> from sympy.physics.continuum_mechanics.beam import Beam
>>> from sympy import symbols
>>> E, I = symbols('E, I')
>>> R1, M1 = symbols('R1, M1')
>>> b = Beam(9, E, I)
>>> b.apply_load(R1, 0, -1)
>>> b.apply_load(M1, 0, -2)
>>> b.apply_load(-8, 0, 0, end=5)
>>> b.apply_load(50, 5, -2)
>>> b.apply_load(-12, 9, -1)
>>> b.bc_slope.append((0, 0))
>>> b.bc_deflection.append((0, 0))
>>> b.solve_for_reaction_loads(R1, M1)
>>> b.reaction_loads
    {M₁: -258, R₁: 52}
>>> b.load
             -2         -1        0             -2            0             -1
    - 258⋅<x>   + 52⋅<x>   - 8⋅<x>  + 50⋅<x - 5>   + 8⋅<x - 5>  - 12⋅<x - 9>
>>> b.shear_force()
             -1         0        1             -1            1             0
    - 258⋅<x>   + 52⋅<x>  - 8⋅<x>  + 50⋅<x - 5>   + 8⋅<x - 5>  - 12⋅<x - 9>
>>> b.bending_moment()
             0         1        2             0            2             1
    - 258⋅<x>  + 52⋅<x>  - 4⋅<x>  + 50⋅<x - 5>  + 4⋅<x - 5>  - 12⋅<x - 9>
>>> b.slope()
                                3                          3             
             1         2   4⋅<x>              1   4⋅<x - 5>             2
    - 258⋅<x>  + 26⋅<x>  - ────── + 50⋅<x - 5>  + ────────── - 6⋅<x - 9> 
                             3                        3                  
    ─────────────────────────────────────────────────────────────────────
                                     E⋅I                                 
>>> b.deflection()
                       3      4                        4             
             2   26⋅<x>    <x>              2   <x - 5>             3
    - 129⋅<x>  + ─────── - ──── + 25⋅<x - 5>  + ──────── - 2⋅<x - 9> 
                    3       3                      3                 
    ─────────────────────────────────────────────────────────────────
                                   E⋅I                               

Example 2
---------

There is a beam of length 30 meters. A moment of magnitude 120 Nm is
applied in the clockwise direction at the end of the beam. A pointload
of magnitude 8 N is applied from the top of the beam at the starting
point. There are two simple supports below the beam. One at the end
and another one at a distance of 10 meters from the start. The
deflection is restricted at both the supports.

::

  || 8 N                                       ⭯ 120 Nm
  \/______________________________________________|
  |_______________________________________________|
              /\                                 /\
  |------------|---------------------------------|
      10 m                  20 m

.. note::

    Using the sign convention of upward forces and clockwise moment
    being positive.

>>> from sympy.physics.continuum_mechanics.beam import Beam
>>> from sympy import symbols
>>> E, I = symbols('E, I')
>>> R1, R2 = symbols('R1, R2')
>>> b = Beam(30, E, I)
>>> b.apply_load(-8, 0, -1)
>>> b.apply_load(R1, 10, -1)
>>> b.apply_load(R2, 30, -1)
>>> b.apply_load(120, 30, -2)
>>> b.bc_deflection.append((10, 0))
>>> b.bc_deflection.append((30, 0))
>>> b.solve_for_reaction_loads(R1, R2)
>>> b.reaction_loads
    {R₁: 6, R₂: 2}
>>> b.load
           -1             -1               -2             -1
    - 8⋅<x>   + 6⋅<x - 10>   + 120⋅<x - 30>   + 2⋅<x - 30>  
>>> b.shear_force()
           0             0               -1             0
    - 8⋅<x>  + 6⋅<x - 10>  + 120⋅<x - 30>   + 2⋅<x - 30> 
>>> b.bending_moment()
           1             1               0             1
    - 8⋅<x>  + 6⋅<x - 10>  + 120⋅<x - 30>  + 2⋅<x - 30> 
>>> b.slope()
           2             2               1           2   4000
    - 4⋅<x>  + 3⋅<x - 10>  + 120⋅<x - 30>  + <x - 30>  + ────
                                                          3  
    ─────────────────────────────────────────────────────────
                               E⋅I                           
>>> b.deflection()
                  3                                      3        
    4000⋅x   4⋅<x>            3              2   <x - 30>         
    ────── - ────── + <x - 10>  + 60⋅<x - 30>  + ───────── - 12000
      3        3                                     3            
    ──────────────────────────────────────────────────────────────
                                 E⋅I                              

Example 3
---------

A beam of length 6 meters is having a roller support at the start and a
hinged support at the end. A clockwise moment of 1.5 kN-m is applied at the mid
of the beam. A constant distributed load of 3 kN/m and a ramp load of 1 kN/m
is applied from the mid till the end of the beam.

::

                              ramp load = 1 KN/m
                            constant load = 3 KN/m
                         |------------------------|
                       ⭯ 1.5 KN-m               
   ______________________|________________________
  |_______________________________________________|
  o                      |                       /\
  |----------------------|-----------------------|
          3.0 m                     3.0 m

.. note::

    Using the sign convention of upward forces and clockwise moment
    being positive.

>>> from sympy.physics.continuum_mechanics.beam import Beam
>>> from sympy import symbols
>>> E, I = symbols('E, I')
>>> R1, R2 = symbols('R1, R2')
>>> b = Beam(6, E, I)
>>> b.apply_load(R1, 0, -1)
>>> b.apply_load(1.5, 3, -2)
>>> b.apply_load(-3, 3, 0)
>>> b.apply_load(-1, 3, 1)
>>> b.apply_load(R2, 6, -1)
>>> b.bc_deflection.append((0, 0))
>>> b.bc_deflection.append((6, 0))
>>> b.solve_for_reaction_loads(R1, R2)
>>> b.reaction_loads
    {R₁: 2.75, R₂: 10.75}
>>> b.load
            -1              -2            0          1                -1
    2.75⋅<x>   + 1.5⋅<x - 3>   - 3⋅<x - 3>  - <x - 3>  + 10.75⋅<x - 6>  
>>> b.shear_force()
                                                    2                 
            0              -1            1   <x - 3>                 0
    2.75⋅<x>  + 1.5⋅<x - 3>   - 3⋅<x - 3>  - ──────── + 10.75⋅<x - 6> 
                                                2                     
>>> b.bending_moment()
                                        2          3                 
            1              0   3⋅<x - 3>    <x - 3>                 1
    2.75⋅<x>  + 1.5⋅<x - 3>  - ────────── - ──────── + 10.75⋅<x - 6> 
                                   2           6                     
>>> b.slope()
                                       3          4                        
             2              1   <x - 3>    <x - 3>                 2       
    1.375⋅<x>  + 1.5⋅<x - 3>  - ──────── - ──────── + 5.375⋅<x - 6>  - 15.6
                                   2          24                           
    ───────────────────────────────────────────────────────────────────────
                                      E⋅I                                  
>>> b.deflection()
                                                              4          5                            
                                   3               2   <x - 3>    <x - 3>                            3
    -15.6⋅x + 0.458333333333333⋅<x>  + 0.75⋅<x - 3>  - ──────── - ──────── + 1.79166666666667⋅<x - 6> 
                                                          8         120                               
    ──────────────────────────────────────────────────────────────────────────────────────────────────
                                                   E⋅I                                                

Example 4
---------

An overhanging beam of length 8 meters is pinned at 1 meter from starting point
and supported by a roller 1 meter before the other end. It is subjected
to a distributed constant load of 10 KN/m from the starting point till
2 meters away from it. Two pointloads of 20KN and 8KN are applied at
5 meters and 7.5 meters away from the starting point respectively.

::

    10 KN/m
  _____________                 20 KN         8 KN
  | | | | | | |                  |             |
  V V V V V V V                  V             V
   _______________________________________________
  |_______________________________________________|
        /\                                  O
  |-----|------|-----------------|----------|--|--|
     1m    1m          3m              2m   .5m .5m

.. note::

    Using the sign convention of upward forces and clockwise moment
    being positive.

.. code:: pycon

    >>> from sympy.physics.continuum_mechanics.beam import Beam
    >>> from sympy import symbols
    >>> E,I,M,V = symbols('E I M V')
    >>> b = Beam(8, E, I)
    >>> E,I,R1,R2 = symbols('E I R1 R2')
    >>> b.apply_load(R1, 1, -1)
    >>> b.apply_load(R2, 7, -1)
    >>> b.apply_load(10, 0, 0, end=2)
    >>> b.apply_load(20, 5, -1)
    >>> b.apply_load(8, 7.5, -1)
    >>> b.solve_for_reaction_loads(R1, R2)
    >>> b.reaction_loads
    {R₁: -26.0, R₂: -22.0}
    >>> b.load
          0               -1             0             -1               -1              -1
    10⋅<x>  - 26.0⋅<x - 1>   - 10⋅<x - 2>  + 20⋅<x - 5>   - 22.0⋅<x - 7>   + 8⋅<x - 7.5>  

    >>> b.shear_force()
          1               0             1             0               0              0
    10⋅<x>  - 26.0⋅<x - 1>  - 10⋅<x - 2>  + 20⋅<x - 5>  - 22.0⋅<x - 7>  + 8⋅<x - 7.5>

    >>> b.bending_moment()
         2               1            2             1               1              1
    5⋅<x>  - 26.0⋅<x - 1>  - 5⋅<x - 2>  + 20⋅<x - 5>  - 22.0⋅<x - 7>  + 8⋅<x - 7.5> 

    >>> b.bc_deflection = [(1, 0), (7, 0)]
    >>> b.slope()
         3                            3                                                                
    5⋅<x>                2   5⋅<x - 2>              2               2              2                   
    ────── - 13.0⋅<x - 1>  - ────────── + 10⋅<x - 5>  - 11.0⋅<x - 7>  + 4⋅<x - 7.5>  + 28.2916666666667
      3                          3                                                                     
    ───────────────────────────────────────────────────────────────────────────────────────────────────
                                                    E⋅I                                                
    >>> b.deflection()
                              4                                        4             3                                          3                   
                         5⋅<x>                            3   5⋅<x - 2>    10⋅<x - 5>                            3   4⋅<x - 7.5>                    
    28.2916666666667⋅x + ────── - 4.33333333333333⋅<x - 1>  - ────────── + ─────────── - 3.66666666666667⋅<x - 7>  + ──────────── - 28.7083333333333
                           12                                     12            3                                         3                         
    ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                          E⋅I

Example 5
---------

A cantilever beam of length 6 meters is under downward distributed constant
load with magnitude of 4.0 KN/m from starting point till 2 meters away
from it. A ramp load of 1 kN/m applied from the mid till the end of
the beam. A point load of 12KN is also applied in same direction 4 meters
away from start.

::

                                       .
                                     . |
                           12 KN   . | |
                             |   . | | |
                             V . | | | |
  \\\\|   4 KN/m             . | | | | |
  \\\\|___________         . 1 KN/m  | |
  \\\\|| | | | | |       . V V V V V V V
  \\\\|V V V V V V     |---------------|
  \\\\|________________________________
  \\\\|________________________________|
  \\\\|          :          :          :
  \\\\|----------|-----|----|----------|
          2.0 m     1m   1m      2.0 m

.. note::

    Using the sign convention of upward forces and clockwise moment
    being positive.

.. code:: pycon

    >>> from sympy.physics.continuum_mechanics.beam import Beam
    >>> from sympy import symbols
    >>> E,I,M,V = symbols('E I M V')
    >>> b = Beam(6, E, I)
    >>> b.apply_load(V, 0, -1)
    >>> b.apply_load(M, 0, -2)
    >>> b.apply_load(4, 0, 0, end=2)
    >>> b.apply_load(12, 4, -1)
    >>> b.apply_load(1, 3, 1, end=6)
    >>> b.solve_for_reaction_loads(V, M)
    >>> b.reaction_loads
    {M: 157/2, V: -49/2}
    >>> b.load
           -2         -1                                                                      
    157⋅<x>     49⋅<x>          0            0          1             -1          0          1
    ───────── - ──────── + 4⋅<x>  - 4⋅<x - 2>  + <x - 3>  + 12⋅<x - 4>   - <x - 6>  - <x - 6> 
        2          2                                                                          
    >>> b.shear_force()
           -1         0                                2                                   2
    157⋅<x>     49⋅<x>         1            1   <x - 3>              0          1   <x - 6> 
    ───────── - ─────── + 4⋅<x>  - 4⋅<x - 2>  + ──────── + 12⋅<x - 4>  - <x - 6>  - ────────
        2          2                               2                                   2    
    >>> b.bending_moment()
           0         1                                3                        2          3
    157⋅<x>    49⋅<x>         2            2   <x - 3>              1   <x - 6>    <x - 6> 
    ──────── - ─────── + 2⋅<x>  - 2⋅<x - 2>  + ──────── + 12⋅<x - 4>  - ──────── - ────────
       2          2                               6                        2          6    
    >>> b.bc_deflection = [(0, 0)]
    >>> b.bc_slope = [(0, 0)]
    >>> b.slope()
           1         2        3            3          4                       3          4
    157⋅<x>    49⋅<x>    2⋅<x>    2⋅<x - 2>    <x - 3>             2   <x - 6>    <x - 6> 
    ──────── - ─────── + ────── - ────────── + ──────── + 6⋅<x - 4>  - ──────── - ────────
       2          4        3          3           24                      6          24   
    ──────────────────────────────────────────────────────────────────────────────────────
                                             E⋅I                                          
    >>> b.deflection()
           2         3      4          4          5                       4          5
    157⋅<x>    49⋅<x>    <x>    <x - 2>    <x - 3>             3   <x - 6>    <x - 6> 
    ──────── - ─────── + ──── - ──────── + ──────── + 2⋅<x - 4>  - ──────── - ────────
       4          12      6        6         120                      24        120   
    ──────────────────────────────────────────────────────────────────────────────────
                                           E⋅I                                        

Example 6
---------

An overhanging beam of length 11 meters is subjected to a distributed constant
load of 2 KN/m from 2 meters away from the starting point till 6 meters away
from it. It is pinned at the starting point and is resting over a roller 8 meters
away from that end. Also a clockwise moment of 5 KN-m is applied at the overhanging end.

::

                 2 KN/m
             _________________
             | | | | | | | | |
             V V V V V V V V V                        ⭯ 5 KN-m
    ____________________________________________________|
   O____________________________________________________|
  / \                                   /\
   |--------|----------------|----------|---------------|
       2m           4m            2m            3m

.. note::

    Using the sign convention of upward forces and clockwise moment
    being positive.

.. code:: pycon

    >>> R1, R2 = symbols('R1, R2')
    >>> E, I = symbols('E, I')
    >>> b = Beam(11, E, I)
    >>> b.apply_load(R1, 0, -1)
    >>> b.apply_load(2, 2, 0, end=6)
    >>> b.apply_load(R2, 8, -1)
    >>> b.apply_load(5, 11, -2)
    >>> b.solve_for_reaction_loads(R1, R2)
    >>> b.reaction_loads
    {R₁: -37/8, R₂: -27/8}
    >>> b.load
            -1                                       -1               
      37⋅<x>              0            0   27⋅<x - 8>               -2
    - ──────── + 2⋅<x - 2>  - 2⋅<x - 6>  - ──────────── + 5⋅<x - 11>  
         8                                      8                     
    >>> b.shear_force()
            0                                       0               
      37⋅<x>             1            1   27⋅<x - 8>              -1
    - ─────── + 2⋅<x - 2>  - 2⋅<x - 6>  - ─────────── + 5⋅<x - 11>  
         8                                     8                    
    >>> b.bending_moment()
            1                                   1              
      37⋅<x>           2          2   27⋅<x - 8>              0
    - ─────── + <x - 2>  - <x - 6>  - ─────────── + 5⋅<x - 11> 
         8                                 8                   
    >>> b.bc_deflection = [(0, 0), (8, 0)]
    >>> b.slope()
            2          3          3             2                   
      37⋅<x>    <x - 2>    <x - 6>    27⋅<x - 8>              1     
    - ─────── + ──────── - ──────── - ─────────── + 5⋅<x - 11>  + 36
         16        3          3            16                       
    ────────────────────────────────────────────────────────────────
                                  E⋅I                               
    >>> b.deflection()
                 3          4          4            3             2
           37⋅<x>    <x - 2>    <x - 6>    9⋅<x - 8>    5⋅<x - 11> 
    36⋅x - ─────── + ──────── - ──────── - ────────── + ───────────
              48        12         12          16            2     
    ───────────────────────────────────────────────────────────────
                                  E⋅I  
