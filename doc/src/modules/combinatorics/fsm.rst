Finite State Machine
====================

.. module:: sympy.combinatorics.rewritingsystem_fsm

Finite State Machine Representation
-----------------------------------

Given:

- Type of state symbols $N$
- Type of input symbols $T$

We make two definitions of finite state machines below:

Nondeterministic Finite State Machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A transition table for nondeterministic finite state machine is defined as a
mapping:

.. math::
    \delta: N \to (T \to 2^N)

Where each state symbols maps to a mapping between input symbols and the set
of state symbols.

If the current state is $A \in N$ and the input alphabet is $a \in T$,
and if $B \in \delta(A)(a)$, then it is possible to make a transition
from $A$ to $B$ in the finite state machine when the input symbol is $a$.

Similarly, a transition table for nondeterministic finite state machine with
empty input ($\epsilon$) is defined as a mapping:

.. math::
    \delta: N \to ((T | \epsilon) \to 2^N)

Where the singleton type for $\epsilon$ is appended to the types of
input symbols.

If the current state is $A \in N$ and if $B \in \delta(A)(\epsilon)$,
it is possible to make a transition from $A$ to $B$ without taking any input
symbol.

Deterministic Finite State Machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A transition table for deterministic finite state machine is defined as a
mapping:

.. math::
    \delta: N \to (T \to N)

If $A \in N$ is the current state symbol and if $a \in T$ is the input symbol,
$a$ should exist in the domain of $\delta(A)$ for it to be possible to make any
transition. And for each $a$ in the domain, there is only one symbol
$B = \delta(A)(a)$ for the transition.

NFA to DFA Conversion
---------------------

.. autofunction:: epsilon_closure

.. autofunction:: instructions

.. autofunction:: subset_construction

.. autofunction:: subset_construction_epsilon

TODO
----

There is an educational framework of finite state machine to have reference
implementation like
http://www.cs.um.edu.mt/gordon.pace/Research/Software/Relic/

If you want to fill in some nontrivial computation of finite state machine to
use in generic symbolic computation, things can be added here.
