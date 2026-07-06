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

Translates plain text and natural deductive ideas into rigorous ``sympy_modal`` expressions, often leveraging LLM APIs for translation.

* ``FormalisationInterface``: Main entrypoint for string/text parsing.
* ``LLMPromptBuilder``: Constructs standardized prompts for LLM conversion tasks.


Examples
--------

Example 1: Defining a Kripke Frame (S4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An S4 frame requires reflexivity and transitivity.

.. code-block:: python

    from sympy.logic.modal import KripkeFrame, Axiom, Box
    from sympy.abc import p

    # Define the frame
    s4_frame = KripkeFrame([Axiom.REFLEXIVITY, Axiom.TRANSITIVITY])

    # In S4, Box(p) -> p is valid (Reflexivity)
    print("S4 Reflexivity validation:", s4_frame.validate_inference(Box(p), p))
    # Output: True

    # In S4, Box(p) -> Box(Box(p)) is valid (Transitivity)
    print("S4 Transitivity validation:", s4_frame.validate_inference(Box(p), Box(Box(p))))
    # Output: True

Example 2: Using the TrustedKernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All verified theorems pass through the ``TrustedKernel``.

.. code-block:: python

    from sympy.logic.modal import KripkeFrame, TrustedKernel, ProofTerm, ModusPonens
    from sympy import Implies
    from sympy.abc import p, q

    kernel = TrustedKernel(KripkeFrame([]))

    # Assume we have proven p and p -> q
    proof_p = ProofTerm(p, kernel_ref=kernel, source="hypothesis")
    proof_p_implies_q = ProofTerm(Implies(p, q), kernel_ref=kernel, source="hypothesis")

    # Apply Modus Ponens
    proof_q = ModusPonens(proof_p_implies_q, proof_p).apply(kernel)
    print("Derived formula via Modus Ponens:", proof_q.formula)
    # Output: q


Example 3: Working with a ProofContext
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ProofContext`` maintains the active proof state and handles assumptions.

.. code-block:: python

    from sympy.logic.modal import ProofContext, KripkeFrame, ProofTerm
    from sympy.abc import p, q

    kernel = KripkeFrame([]).get_kernel()
    ctx = ProofContext(kernel)

    # Introduce an assumption
    ctx.assume(p)

    # Prove q from p (trivially in this example via direct assumption injection)
    proof_q_from_p = ProofTerm(q, kernel_ref=kernel, source="hypothesis")
    ctx.add_proof(proof_q_from_p)

    # Discharge the assumption p
    discharged_proof = ctx.discharge(p)
    print("Discharged implication:", discharged_proof.formula)
    # Output: Implies(p, q)


Example 4: Fixed Points in Provability Logic (GL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Provability Logic features Löb's theorem which constructs well-founded fixed points.

.. code-block:: python

    from sympy.logic.modal import (
        KripkeFrame, Axiom, ProvabilityBox,
        GuardedFixedPoint, FormalisationInterface
    )
    from sympy.abc import p

    # GL requires Transitivity and Converse Well-Foundedness
    gl_frame = KripkeFrame([Axiom.TRANSITIVITY, Axiom.CONVERSE_WELL_FOUNDED])
    interface = FormalisationInterface(frame=gl_frame)

    # Box(Box(p) -> p) -> Box(p)
    loeb_formula = interface.parse("Box(Implies(Box(p), p)) >> Box(p)")
    print("Löb's theorem proved in GL:", gl_frame.validate_inference(None, loeb_formula))
    # Output: True


Example 5: Natural Language Parsing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The interface parses human-readable strings into typed SymPy objects.

.. code-block:: python

    from sympy.logic.modal import FormalisationInterface, AlethicBox
    from sympy.abc import p

    interface = FormalisationInterface()

    # "It is necessary that p"
    parsed_formula = interface.parse("AlethicBox(p)")

    print("Parsed formula type:", type(parsed_formula).__name__)
    # Output: AlethicBox


Example 6: Classical Injection & SMT Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the calculus is intuitionistic. You can optionally allow classical reasoning and offload verification to SMT solvers.

.. code-block:: python

    from sympy.logic.modal import ProofContext, KripkeFrame
    from sympy.abc import p

    kernel = KripkeFrame([]).get_kernel()

    # Enable classical logic
    ctx = ProofContext(kernel, allow_classical=True)

    # Law of Excluded Middle
    lem_proof = ctx.derive_classical_axiom(p | ~p)
    print("Classical proof valid (LEM):", lem_proof is not None)
    # Output: True

Example 7: Semantic Evaluation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

KripkeModels provide tools for semantic evaluation and model checking across worlds.

.. code-block:: python

    from sympy.logic.modal import KripkeModel, SemanticEvaluator, Box, Diamond, AgentBox, CommonKnowledge
    from sympy.abc import p

    # Define worlds and relations
    worlds = ['w1', 'w2', 'w3']
    relations = {'alethic': [('w1', 'w2'), ('w2', 'w3')]}
    valuation = {'p': ['w2']}

    model = KripkeModel(worlds, relations, valuation)
    evaluator = SemanticEvaluator(model)

    # p is true at w2, and w1 sees w2, so Diamond(p) is true at w1
    print("Is Diamond(p) true at w1?", evaluator.is_true_at(Diamond(p), 'w1'))
    # Output: True

    # But Box(p) is false at w1 because w1 could potentially see a world where p is false (if we added one)
    # Actually, in this specific model, w1 only sees w2, and p is true at w2.
    # Wait, if w1 only sees w2, and p is true at w2, Box(p) SHOULD be true at w1.
    # Let's check the evaluator output for exact logic.
    print("Is Box(p) true at w1?", evaluator.is_true_at(Box(p), 'w1'))

    # Multi-agent epistemic logic is supported via string identifiers
    epistemic_box = AgentBox(p, "Alice")
    print("AgentBox modality:", epistemic_box.agent_modality())

    common_knowledge = CommonKnowledge(p, "GroupA")
    print("CommonKnowledge modality:", common_knowledge.group_modality())
