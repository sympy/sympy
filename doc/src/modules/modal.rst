SymPy Modal
===========

``sympy.logic.modal`` is a second-order modal predicate calculus extension to SymPy. It provides a research-grade system combining a computer algebra interface with a proof-theoretic kernel.

Reference Documentation
-----------------------

Core Types & Universe
^^^^^^^^^^^^^^^^^^^^^

The type system replaces SymPy's untyped predicate logic with a cumulative type hierarchy.

* ``Universe(level: int)``: Defines a universe level. ``Universe(0)`` is predicative.
* ``BoolType()``: The type of standard Boolean propositions.
* ``FunctionType(domain, codomain)``: Represents functions (e.g. predicates are functions returning ``BoolType()``).
* ``PredicateVariable(name, type)``: A second-order variable ranging over given types.
* ``ModalPredicate(name, type, args)``: A typed predicate applied to arguments whose validity depends on the frame.
* ``GuardedFixedPoint(operator, name)``: Well-founded fixed-point constructors for modal arguments (Löb's theorem).

Logical Operators
^^^^^^^^^^^^^^^^^

``sympy.logic.modal`` provides the following logical primitives:

* ``Box(formula)``: Necessity (□).
* ``Diamond(formula)``: Possibility (◇). Equivalent to ``~Box(~P)``.
* ``ForAllPredicates(variable, formula)``: Second-order universal quantification.
* ``ExistsPredicates(variable, formula)``: Second-order existential quantification.

For parsing ease, there are named modalities that imply their logic frame automatically via the ``FormalisationInterface``:

* ``ProvabilityBox`` (GL)
* ``AlethicBox`` (S5)
* ``EpistemicBox`` (K45)
* ``DeonticBox`` (D)
* ``TemporalBox`` (S4)

Kripke Frame Semantics
^^^^^^^^^^^^^^^^^^^^^^

Kripke frames define the properties of the accessibility relation used by the modal operators.

* ``KripkeFrame(axioms)``: Instantiates a frame constrained by the given axioms.
* ``Axiom``: An Enum detailing known properties (e.g., ``Axiom.REFLEXIVITY``, ``Axiom.TRANSITIVITY``).

Proof-Theoretic Kernel
^^^^^^^^^^^^^^^^^^^^^^

The ``TrustedKernel`` is the TCB (Trusted Computing Base) for inferences.

* ``TrustedKernel(frame)``: Verifies inference rules and creates proof objects based on frame axioms.
* ``ModusPonens(premise1, premise2)``: Reifies the modus ponens rule.
* ``ProofContext(kernel)``: Manages assumptions, derived proofs, and handles proof search via ``Strategy`` algorithms.

Formalisation Bridge
^^^^^^^^^^^^^^^^^^^^

Translates plain text and natural deductive ideas into rigorous ``sympy.logic.modal`` expressions.
.. , often leveraging LLM APIs for translation.

* ``FormalisationInterface``: Main entrypoint for string/text parsing.

.. * ``LLMPromptBuilder``: Constructs standardized prompts for LLM conversion tasks.

.. note::
   The ``LLMPromptBuilder`` functionality for translating natural language to ``sympy.logic.modal``
   expressions using OpenRouter and Gemini 2.5 Flash has been intentionally commented out of the
   source code and documentation.

   If you wish to experiment with this feature in the future, you can uncomment it in the following files:

   * ``sympy/logic/modal/formalise.py`` (class definition and API requests)
   * ``sympy/logic/modal/__init__.py`` (module exports)
   * ``sympy/__init__.py`` (top-level exports)
   * ``sympy/logic/modal/tests/test_formalise_bridge.py`` (tests)
   * ``doc/src/modules/modal.rst`` (this documentation)


Examples
--------

Example 1: Defining a Kripke Frame (S4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An S4 frame requires reflexivity and transitivity.

.. code-block:: python

    from sympy import Symbol, Implies
    from sympy.logic.modal import KripkeFrame, Box

    # Define the S4 frame
    s4 = KripkeFrame.S4()
    p = Symbol('p')

    # In S4, Box(p) -> p is valid (Reflexivity)
    print("S4 Reflexivity validation:", s4.validates(Implies(Box(p), p)))
    # Output: True

    # In S4, Box(p) -> Box(Box(p)) is valid (Transitivity)
    print("S4 Transitivity validation:", s4.validates(Implies(Box(p), Box(Box(p)))))
    # Output: True

Example 2: Using the TrustedKernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All verified theorems pass through the ``TrustedKernel``.

.. code-block:: python

    from sympy import Symbol, Implies
    from sympy.logic.modal import KripkeFrame, TrustedKernel, ProofTerm, ModusPonens

    kernel = TrustedKernel(frame=KripkeFrame.GL())
    p = Symbol('p')
    q = Symbol('q')

    # Assume we have proven p and p -> q
    pt_ab = ProofTerm(Implies(p, q, evaluate=False), source="axiom")
    pt_a = ProofTerm(p)

    # Apply Modus Ponens verification
    pt_b = kernel.verify_rule(ModusPonens, [pt_ab, pt_a])
    print("Derived formula via Modus Ponens:", pt_b.formula)
    # Output: q


Example 3: Working with a ProofContext
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ProofContext`` maintains the active proof state and handles assumptions.

.. code-block:: python

    from sympy import Symbol, Implies
    from sympy.logic.modal import ProofContext, KripkeFrame, ProofTerm

    ctx = ProofContext(KripkeFrame.K())
    p = Symbol('p')
    q = Symbol('q')

    # Introduce an assumption
    hyp = ctx.assume(p)
    print("Hypothesis source:", hyp.source)

    # Create a fake derivation from the hypothesis for demonstration
    pt_q = ProofTerm(q, hypotheses=[p])

    # Discharging the hypothesis yields P -> Q
    impl = ctx.discharge(hyp, pt_q)
    print("Discharged implication:", impl.formula)
    # Output: Implies(p, q)


Example 4: Fixed Points in Provability Logic (GL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Provability Logic features Löb's theorem which constructs well-founded fixed points.

.. code-block:: python

    from sympy import symbols, Implies
    from sympy.logic.modal import (
        ProofContext, KripkeFrame, PredicateVariable,
        ForAllPredicates, Box, Universe, FunctionType, BoolType
    )

    x = symbols('x')
    ctx = ProofContext(frame=KripkeFrame.GL())
    P = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))

    lob = ForAllPredicates(P,
        Implies(Box(Implies(Box(P(x)), P(x))), Box(P(x)))
    )

    # Certificate proves the theorem is valid under GL frame
    proof = ctx.prove(lob)
    print("Löb's theorem proved in GL:", proof.is_valid)
    # Output: True


Example 5: Natural Language Parsing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The interface parses human-readable strings into typed SymPy objects.

.. code-block:: python

    from sympy.logic.modal import FormalisationInterface, AlethicBox

    fi = FormalisationInterface()

    # "It is necessary that p"
    expr_str = "AlethicBox(Symbol('p'))"
    formula = fi.formalise(expr_str)

    print("Parsed formula type:", type(formula).__name__)
    # Output: AlethicBox

    # The interface can infer the correct frame automatically
    sig = fi.resolve_modality(expr_str)
    print("Inferred modality signatures:", sig.operators)


Example 6: Classical Injection & SMT Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the calculus is intuitionistic. You can optionally allow classical reasoning and offload verification to SMT solvers.

.. code-block:: python

    from sympy import Symbol, Or, Not, Implies
    from sympy.logic.modal import ProofContext, KripkeFrame

    # Standard context rejects LEM without the flag
    ctx_intuitionistic = ProofContext(KripkeFrame.K())
    p = Symbol('p')
    lem = Or(p, Not(p))
    proof_fail = ctx_intuitionistic.prove(lem)
    print("Intuitionistic proof valid:", isinstance(proof_fail, type(ctx_intuitionistic.assume(p))))
    # Output: False

    # With classical logic enabled
    ctx_classical = ProofContext(KripkeFrame.K(), allow_classical=True)
    proof_success = ctx_classical.prove(lem)
    print("Classical proof valid (LEM):", proof_success.is_valid)
    print("Classical derivation source:", proof_success.derivation)
    # Output: True, ('SMT_Solver', p | ~p)

    # SMT solver solving a complex propositional tautology instantly
    q = Symbol('q')
    complex_tautology = Or(Implies(p, q), Implies(q, p))
    smt_proof = ctx_classical.prove(complex_tautology)
    print("SMT proof valid:", smt_proof.is_valid)


Example 7: Semantic Evaluation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

KripkeModels provide tools for semantic evaluation and model checking across worlds.

.. code-block:: python

    from sympy import Symbol
    from sympy.logic.modal import KripkeModel, SemanticEvaluator, Box, Diamond, AgentBox, CommonKnowledge

    # Kripke Model setup
    W = {'w1', 'w2', 'w3'}
    R = {'w1': {'w2', 'w3'}}
    p = Symbol('p')
    V = {
        ('w2', p): True,
        ('w3', p): False
    }
    model = KripkeModel(W, R, V)
    evaluator = SemanticEvaluator(model)

    # Evaluation
    print("Is Diamond(p) true at w1?", evaluator.evaluate(Diamond(p), 'w1'))
    # Output: True

    print("Is Box(p) true at w1?", evaluator.evaluate(Box(p), 'w1'))
    # Output: False

    # Expressive Multi-Agent Operators
    agent_box = AgentBox('Alice', p)
    common_knowledge = CommonKnowledge('GroupA', p)
    print("AgentBox modality:", agent_box.modality)
    print("CommonKnowledge modality:", common_knowledge.modality)
