# Analysis of Algebraic Tensor Implementation Structure

## Executive Summary

This document analyzes the current class hierarchy design in the algebraic tensor implementation and provides recommendations for restructuring to align with SymPy's core architecture for successful mainline integration.

### Current Implementation

```
AlgebraicPureTensor(Mul)
AlgebraicTensor(Basic)
AlgebraicZeroTensor(AtomicExpr)
```

### Recommended Structure

```
All three classes should extend Basic directly
```

---

## 1. SymPy Core Class Hierarchy Reference

### 1.1 Basic -> Expr -> AtomicExpr

```
Basic (core/basic.py)
    |
    +-- Atom (atomic objects with no subexpressions)
    |       |
    |       +-- AtomicExpr (atoms that are also expressions)
    |
    +-- Expr (algebraic expressions with arithmetic operations)
```

**AtomicExpr Key Properties:**
- `is_Atom = True` - Marks atomic nature (no `.args` to traverse)
- `_eval_derivative(self, s)` - Returns `S.One` if `self == s`, else `S.Zero`
- Designed for indivisible symbolic objects (Numbers, Symbols, etc.)
- No internal structure for pattern matching or simplification

### 1.2 Basic -> AssocOp -> Mul

```
Basic (core/basic.py)
    |
    +-- AssocOp (core/operations.py) - associative operations base class
            |
            +-- Mul (core/mul.py) - multiplication with flatten()
            +-- Add (core/add.py) - addition with similar structure
```

**AssocOp Key Properties:**
- `__slots__ = ('is_commutative',)` - Commutativity stored directly for performance
- `_from_args(cls, args, is_commutative=None)` - Creates instance from pre-processed args
- `flatten(cls, seq)` - Abstract method for argument normalization
- Handles associativity, identity elimination, single-arg unwrapping

**Mul Specific Design:**
- `is_Mul = True` - Type flag
- `_args_type = Expr` - Enforces Expr arguments
- `flatten()` is the heart of evaluation: separates commutative/non-commutative parts
- `is_commutative` computed from presence of non-commutative args
- Coefficient stored as first argument

---

## 2. Analysis of Current Design Decisions

### 2.1 AlgebraicPureTensor extending Mul

**Current Rationale (from code comments):**
> "Extends SymPy's non-commutative Mul so that all existing simplification, differentiation, and pattern-matching machinery is available out of the box."

**Problems:**

1. **`flatten()` Incompatibility**
   - `Mul.flatten()` expects to separate commutative coefficients from non-commutative factors using standard SymPy rules
   - `AlgebraicPureTensor` stores coefficients as first arg with factors following, but `Mul.flatten()` will attempt to:
     - Combine numeric powers (`x**a * x**b -> x**(a+b)`)
     - Apply numeric simplification (`2*3 -> 6`)
     - Order factors canonically
   - These transformations may break tensor-specific semantics

2. **`is_commutative` Conflict**
   - `Mul` computes `is_commutative = not nc_part` from flattened args
   - `AlgebraicPureTensor` overrides `_eval_is_commutative = lambda self: False` unconditionally
   - This creates inconsistency: the base class computes one value, the subclass asserts another
   - The `is_Mul = False` hack to prevent `Mul.flatten()` from unpacking shows awareness of this conflict

3. **Type System Pollution**
   - `isinstance(pt, Mul)` returns `True`, which may confuse code expecting standard multiplication
   - `Mul`'s `args_cnc()`, `as_coeff_Mul()`, and other methods may not work correctly with tensor factors

4. **Evaluation Control Issues**
   - `Mul.__new__` calls `flatten()` unconditionally unless `evaluate=False`
   - The `flatten()` method is complex and may introduce unwanted transformations
   - The current implementation works around this with careful `__new__` logic

**What Actually Works:**
- The `__slots__ = ()` inheritance is compatible
- `__new__` bypasses `flatten()` by calling `Mul.__new__(cls, ...)` directly with `evaluate=False`
- The `is_Mul = False` flag prevents some unwanted behavior

**Better Alternative:**
- Extend `Basic` directly and implement tensor-specific composition logic
- Reuse `Mul` only for coefficient arithmetic, not as a base class

---

### 2.2 AlgebraicTensor extending Basic

**Current Rationale:**
> "Built on top of Basic (not Add) to avoid the is_commutative descriptor conflict in AssocOp._from_args, while still exposing is_Add = True so that the wider SymPy ecosystem recognises this as an additive expression."

**Analysis:**

1. **Correct Assessment**
   - The comment correctly identifies the `is_commutative` descriptor conflict
   - `AssocOp._from_args()` sets `obj.is_commutative` as a slot attribute
   - If `AlgebraicTensor` extended `Add`, this would conflict with tensor-specific commutativity logic

2. **`is_Add = True` Hack**
   - Setting `is_Add = True` signals to SymPy that this behaves like addition
   - However, this is incomplete: `AlgebraicTensor` does not inherit Add's simplification, canonicalization, or pattern-matching
   - External code expecting `Add` semantics may not work correctly

3. **Missing Add Machinery**
   - No `flatten()` method for term collection and coefficient combination
   - No canonical ordering of terms
   - No automatic simplification (`x + (-x) -> 0`)
   - Manual implementation of `_merge()`, `_subtract()`, `_collect_coefficients()` duplicates Add logic

**Problems:**

1. **Inconsistent with SymPy Conventions**
   - Sums of terms typically extend `Add` or at least follow `Add` patterns
   - Pattern matching against `Add` may not match `AlgebraicTensor`
   - `as_coeff_Add()`, `as_coefficients_dict()` may not work

2. **Maintenance Burden**
   - Re-implementing Add-like functionality creates divergence from SymPy's evolution
   - Bug fixes and improvements to `Add` won't propagate

**Better Alternative:**
- Extending `Basic` is the correct choice, but consider whether `AlgebraicTensor` needs the full Add machinery
- If tensor sums have unique semantics, `Basic` is appropriate
- Document the divergence clearly and provide compatibility adapters

---

### 2.3 AlgebraicZeroTensor extending AtomicExpr

**Current Rationale:**
> "Extends AtomicExpr (Atom + Expr) to integrate with SymPy's expression system: sympify(), tree traversal (atoms(), has(), replace()), the assumptions system (is_commutative), generic operations (subs(), xreplace(), doit()), and the Expr machinery (as_base_exp(), as_coeff_Mul(), as_coeff_Add())."

**Critical Problems:**

1. **Fundamental Misunderstanding of AtomicExpr**
   - `AtomicExpr` is for **atoms** - objects with **no internal structure**
   - `AtomicExpr` inherits from `Atom`, which has `is_Atom = True`
   - Atoms are indivisible: `atoms()` returns the atom itself, not subcomponents
   - Pattern matching treats atoms as leaves

2. **AlgebraicZeroTensor Has Structure**
   - It has a `_shape` attribute that defines its identity
   - Different shapes produce different zero tensors
   - The shape is integral to the object's semantics
   - This is **not** atomic - it's a structured object

3. **Hash and Equality Issues**
   - `AtomicExpr` uses Basic's default hashing based on type and args
   - `AlgebraicZeroTensor` stores shape in `_shape` (not in `args`)
   - `_hashable_content()` returns `(_shape,)`, but this may conflict with Basic's expectations
   - Objects with different internal storage but same type may hash incorrectly

4. **args Property Semantics**
   - `Basic.args` should return the structural components
   - `AlgebraicZeroTensor.args` returns `()` (empty tuple from AtomicExpr)
   - This misleads code that inspects `.args` for structure

5. **Differentiation Semantics**
   - `AtomicExpr._eval_derivative()` returns `S.One` if `self == s`, else `S.Zero`
   - For `AlgebraicZeroTensor`, the derivative should always be a zero tensor of the same shape
   - The current `diff()` override works, but the AtomicExpr base class semantics are wrong

6. **Assumption System Conflicts**
   - `AtomicExpr` is designed for indivisible symbolic entities
   - `AlgebraicZeroTensor` is a structured zero with shape-dependent behavior
   - The assumptions system may treat it as a leaf when it has internal structure

**What Actually Works:**
- `is_commutative = True` is correctly set
- `is_zero = True` integrates with assumptions
- Manual overrides (`diff()`, `expand()`, `doit()`) provide correct behavior
- The `_shape` attribute works for shape tracking

**Better Alternative:**
- Extend `Basic` directly
- Store shape in `args` or ensure `_hashable_content()` properly includes it
- This makes the structural nature explicit and consistent with SymPy conventions

---

## 3. Recommended High-Level Structure

### 3.1 Unified Basic-Based Hierarchy

```
All classes extend Basic directly:

class AlgebraicPureTensor(Basic):
    # Non-commutative tensor product of factors
    # Extends Basic, not Mul
    
class AlgebraicTensor(Basic):
    # Sum of AlgebraicPureTensors with same shape
    # Already correct, keep as Basic
    
class AlgebraicZeroTensor(Basic):
    # Zero tensor with shape
    # Change from AtomicExpr to Basic
```

### 3.2 Rationale

1. **Consistency**
   - All tensor classes follow the same inheritance pattern
   - No confusing mixins or special-case base classes
   - Clear separation from SymPy's arithmetic operators

2. **Avoid Mul Pitfalls**
   - No unwanted `flatten()` transformations
   - No `is_commutative` descriptor conflicts
   - Full control over evaluation and simplification

3. **Proper Atomic vs. Structural Semantics**
   - `AlgebraicZeroTensor` is not atomic - it has shape structure
   - `Basic` correctly represents structured symbolic objects
   - No confusion about `.args` or tree traversal

4. **Easier Integration**
   - SymPy maintainers can review a clean, self-contained design
   - No need to justify deviations from Mul/Add patterns
   - Clear documentation of tensor-specific semantics

### 3.3 Implementation Guidelines

#### For AlgebraicPureTensor:

```python
class AlgebraicPureTensor(Basic):
    """Pure tensor as a non-commutative tensor product."""
    
    __slots__ = ()
    
    is_AlgebraicPureTensor = True
    # Remove: is_Mul = False (no longer relevant)
    
    _eval_is_commutative = lambda self: False
    
    @classmethod
    def __new__(cls, *args, evaluate=True):
        # Validate arguments
        # Separate coefficient from factors
        # Handle zero cases
        # Call Basic.__new__(cls, *processed_args)
        pass
    
    @property
    def args(self):
        # Return (coeff, *factors) or just *factors
        # Ensure args captures full structure
        pass
```

**Key Changes:**
- Call `Basic.__new__()` instead of `Mul.__new__()`
- Implement custom evaluation logic
- Store coefficient and factors explicitly in args
- No reliance on Mul's flatten()

#### For AlgebraicZeroTensor:

```python
class AlgebraicZeroTensor(Basic):
    """Zero tensor carrying a specific tensor shape."""
    
    __slots__ = ('_shape',)
    
    is_AlgebraicZeroTensor = True
    is_zero = True
    is_commutative = True  # This works with Basic
    
    def __new__(cls, shape):
        # Normalize shape
        # Create instance with Basic.__new__
        # Store shape in _shape slot
        obj = Basic.__new__(cls)
        obj._shape = shape
        return obj
    
    def _hashable_content(self):
        return (self._shape,)
    
    @property
    def shape(self):
        return self._shape
```

**Key Changes:**
- Extend `Basic` instead of `AtomicExpr`
- Keep `_shape` in slot for efficiency
- `_hashable_content()` already correct
- No change needed to most methods

#### For AlgebraicTensor:

```python
class AlgebraicTensor(Basic):
    """Sum of AlgebraicPureTensors sharing the same shape."""
    
    # Current implementation is already correct
    # Keep as Basic extension
    pass
```

### 3.4 Compatibility Considerations

1. **Pattern Matching**
   - Update pattern matching docs to reference `AlgebraicPureTensor`/`AlgebraicTensor` explicitly
   - Not `Mul`/`Add`

2. **Type Checking**
   - Replace `isinstance(x, Mul)` with `isinstance(x, AlgebraicPureTensor)`
   - Replace `isinstance(x, Add)` with `isinstance(x, AlgebraicTensor)`

3. **Sympification**
   - Ensure `converter` dictionary or `_sympy_` methods handle tensor types
   - Test with `sympify()` on tensor expressions

4. **Printing**
   - Verify `_repr_latex_()` and printers work correctly
   - No changes expected (printing is based on `func` and `args`)

---

## 4. Migration Checklist

### Phase 1: AlgebraicZeroTensor → Basic

- [ ] Change base class from `AtomicExpr` to `Basic`
- [ ] Verify `_hashable_content()` and equality
- [ ] Test `sympify()`, `atoms()`, `has()`
- [ ] Test differentiation, substitution
- [ ] Update imports in all modules

### Phase 2: AlgebraicPureTensor → Basic

- [ ] Change base class from `Mul` to `Basic`
- [ ] Rewrite `__new__()` to call `Basic.__new__()`
- [ ] Implement explicit arg structure (coeff, factors)
- [ ] Migrate `__mul__`, `__rmul__`, `__add__`, etc.
- [ ] Test coefficient handling
- [ ] Test `diff()`, `doit()`, `expand()`, `simplify()`
- [ ] Update all internal references

### Phase 3: Integration Testing

- [ ] Run full test suite
- [ ] Test with external SymPy code that uses tensors
- [ ] Verify pattern matching works
- [ ] Verify printing works
- [ ] Test pickling/unpickling

### Phase 4: Documentation

- [ ] Update module docstrings
- [ ] Document inheritance decisions
- [ ] Add migration guide for users
- [ ] Update API reference

---

## 5. Conclusion

The current design shows thoughtful consideration of SymPy's architecture, particularly the decision to avoid `Add` for `AlgebraicTensor` due to the `is_commutative` conflict. However, extending `Mul` for `AlgebraicPureTensor` and `AtomicExpr` for `AlgebraicZeroTensor` introduces subtle incompatibilities that could hinder mainline integration.

**The path forward is clear:**
1. Extend `Basic` for all three classes
2. Implement tensor-specific evaluation logic explicitly
3. Avoid inheriting from arithmetic operators that have strong semantic expectations

This approach provides:
- Clean separation from SymPy's core arithmetic
- Full control over tensor semantics
- Easier maintenance and review
- Better alignment with SymPy's extension patterns

The restructuring effort is moderate but necessary for long-term maintainability and successful integration into the main SymPy codebase.
