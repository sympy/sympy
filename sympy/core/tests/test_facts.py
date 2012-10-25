from sympy.core.facts import (deduce_alpha_implications,
        apply_beta_to_alpha_route, rules_2prereq, FactRules, FactKB)
from sympy.core.logic import And, Not
from sympy.utilities.pytest import XFAIL, raises

T = True
F = False
U = None


def test_deduce_alpha_implications():
    def D(i):
        I = deduce_alpha_implications(i)
        P = rules_2prereq(dict(
            ((k, True), set([(v, True) for v in S])) for k, S in I.items()))
        return I, P

    # transitivity
    I, P = D([('a', 'b'), ('b', 'c')])
    assert I == {'a': set(['b', 'c']), 'b': set(['c'])}
    assert P == {'b': set(['a']), 'c': set(['a', 'b'])}

    # see if the output does not contain repeated implications
    I, P = D([('a', 'b'), ('b', 'c'), ('b', 'c')])
    assert I == {'a': set(['b', 'c']), 'b': set(['c'])}
    assert P == {'b': set(['a']), 'c': set(['a', 'b'])}

    # see if it is tolerant to cycles
    assert D([('a', 'a'), ('a', 'a')]) == ({}, {})
    assert D([('a', 'b'), ('b', 'a')]) == (
        {'a': set(['b']), 'b': set(['a'])},
        {'a': set(['b']), 'b': set(['a'])})

    # see if it catches inconsistency
    raises(ValueError, lambda: D([('a', Not('a'))]))
    raises(ValueError, lambda: D([('a', 'b'), ('b', Not('a'))]))
    raises(ValueError, lambda: D([('a', 'b'), ('b', 'c'), ('b', 'na'),
           ('na', Not('a'))]))

    # something related to real-world
    I, P = D([('rat', 'real'), ('int', 'rat')])

    assert I == {'int': set(['rat', 'real']), 'rat': set(['real'])}
    assert P == {'rat': set(['int']), 'real': set(['int', 'rat'])}


# TODO move me to appropriate place
def test_apply_beta_to_alpha_route():
    APPLY = apply_beta_to_alpha_route

    # indicates empty alpha-chain with attached beta-rule #bidx
    def Q(bidx):
        return (set(), [bidx])

    # x -> a        &(a,b) -> x     --  x -> a
    A = {'x': set(['a'])}
    B = [(And('a', 'b'), 'x')]
    assert APPLY(A, B) == {'x': (set(['a']), []), 'a': Q(0), 'b': Q(0)}

    # x -> a        &(a,!x) -> b    --  x -> a
    A = {'x': set(['a'])}
    B = [(And('a', Not('x')), 'b')]
    assert APPLY(A, B) == {'x': (set(['a']), []), Not('x'): Q(0), 'a': Q(0)}

    # x -> a b      &(a,b) -> c     --  x -> a b c
    A = {'x': set(['a', 'b'])}
    B = [(And('a', 'b'), 'c')]
    assert APPLY(A, B) == \
        {'x': (set(['a', 'b', 'c']), []), 'a': Q(0), 'b': Q(0)}

    # x -> a        &(a,b) -> y     --  x -> a [#0]
    A = {'x': set(['a'])}
    B = [(And('a', 'b'), 'y')]
    assert APPLY(A, B) == {'x': (set(['a']), [0]), 'a': Q(0), 'b': Q(0)}

    # x -> a b c    &(a,b) -> c     --  x -> a b c
    A = {'x': set(['a', 'b', 'c'])}
    B = [(And('a', 'b'), 'c')]
    assert APPLY(A, B) == \
        {'x': (set(['a', 'b', 'c']), []), 'a': Q(0), 'b': Q(0)}

    # x -> a b      &(a,b,c) -> y   --  x -> a b [#0]
    A = {'x': set(['a', 'b'])}
    B = [(And('a', 'b', 'c'), 'y')]
    assert APPLY(A, B) == \
        {'x': (set(['a', 'b']), [0]), 'a': Q(0), 'b': Q(0), 'c': Q(0)}

    # x -> a b      &(a,b) -> c     --  x -> a b c d
    # c -> d                            c -> d
    A = {'x': set(['a', 'b']), 'c': set(['d'])}
    B = [(And('a', 'b'), 'c')]
    assert APPLY(A, B) == {'x': (set(['a', 'b', 'c', 'd']), []),
        'c': (set(['d']), []), 'a': Q(0), 'b': Q(0)}

    # x -> a b      &(a,b) -> c     --  x -> a b c d e
    # c -> d        &(c,d) -> e         c -> d e
    A = {'x': set(['a', 'b']), 'c': set(['d'])}
    B = [(And('a', 'b'), 'c'), (And('c', 'd'), 'e')]
    assert APPLY(A, B) == {'x': (set(['a', 'b', 'c', 'd', 'e']), []),
        'c': (set(['d', 'e']), []), 'a': Q(0), 'b': Q(0), 'd': Q(1)}

    # x -> a b      &(a,y) -> z     --  x -> a b y z
    #               &(a,b) -> y
    A = {'x': set(['a', 'b'])}
    B = [(And('a', 'y'), 'z'), (And('a', 'b'), 'y')]
    assert APPLY(A, B) == {'x': (set(['a', 'b', 'y', 'z']), []),
        'a': (set(), [0, 1]), 'y': Q(0), 'b': Q(1)}

    # x -> a b      &(a,!b) -> c    --  x -> a b
    A = {'x': set(['a', 'b'])}
    B = [(And('a', Not('b')), 'c')]
    assert APPLY(A, B) == \
        {'x': (set(['a', 'b']), []), 'a': Q(0), Not('b'): Q(0)}

    # !x -> !a !b   &(!a,b) -> c    --  !x -> !a !b
    A = {Not('x'): set([Not('a'), Not('b')])}
    B = [(And(Not('a'), 'b'), 'c')]
    assert APPLY(A, B) == \
        {Not('x'): (set([Not('a'), Not('b')]), []), Not('a'): Q(0), 'b': Q(0)}

    # x -> a b      &(b,c) -> !a    --  x -> a b
    A = {'x': set(['a', 'b'])}
    B = [(And('b', 'c'), Not('a'))]
    assert APPLY(A, B) == {'x': (set(['a', 'b']), []), 'b': Q(0), 'c': Q(0)}

    # x -> a b      &(a, b) -> c    --  x -> a b c p
    # c -> p a
    A = {'x': set(['a', 'b']), 'c': set(['p', 'a'])}
    B = [(And('a', 'b'), 'c')]
    assert APPLY(A, B) == {'x': (set(['a', 'b', 'c', 'p']), []),
        'c': (set(['p', 'a']), []), 'a': Q(0), 'b': Q(0)}


def test_FactRules_parse():
    f = FactRules('a -> b')
    assert f.prereq == {'b': set(['a']), 'a': set(['b'])}

    f = FactRules('a -> !b')
    assert f.prereq == {'b': set(['a']), 'a': set(['b'])}

    f = FactRules('!a -> b')
    assert f.prereq == {'b': set(['a']), 'a': set(['b'])}

    f = FactRules('!a -> !b')
    assert f.prereq == {'b': set(['a']), 'a': set(['b'])}

    f = FactRules('!z == nz')
    assert f.prereq == {'z': set(['nz']), 'nz': set(['z'])}


def test_FactRules_parse2():
    raises(ValueError, lambda: FactRules('a -> !a'))


def test_FactRules_deduce():
    f = FactRules(['a -> b', 'b -> c', 'b -> d', 'c -> e'])

    def D(facts):
        kb = FactKB(f)
        kb.deduce_all_facts(facts)
        return kb

    assert D({'a': T}) == {'a': T, 'b': T, 'c': T, 'd': T, 'e': T}
    assert D({'b': T}) == {        'b': T, 'c': T, 'd': T, 'e': T}
    assert D({'c': T}) == {                'c': T,         'e': T}
    assert D({'d': T}) == {                        'd': T        }
    assert D({'e': T}) == {                                'e': T}

    assert D({'a': F}) == {'a': F                                }
    assert D({'b': F}) == {'a': F, 'b': F                        }
    assert D({'c': F}) == {'a': F, 'b': F, 'c': F                }
    assert D({'d': F}) == {'a': F, 'b': F,         'd': F        }

    assert D({'a': U}) == {'a': U}  # XXX ok?


def test_FactRules_deduce2():
    # pos/neg/zero, but the rules are not sufficient to derive all relations
    f = FactRules(['pos -> !neg', 'pos -> !z'])

    def D(facts):
        kb = FactKB(f)
        kb.deduce_all_facts(facts)
        return kb

    assert D({'pos': T}) == {'pos': T, 'neg': F, 'z': F}
    assert D({'pos': F}) == {'pos': F                  }
    assert D({'neg': T}) == {'pos': F, 'neg': T        }
    assert D({'neg': F}) == {          'neg': F        }
    assert D({'z': T}) == {'pos': F,           'z': T}
    assert D({'z': F}) == {                    'z': F}

    # pos/neg/zero. rules are sufficient to derive all relations
    f = FactRules(['pos -> !neg', 'neg -> !pos', 'pos -> !z', 'neg -> !z'])

    assert D({'pos': T}) == {'pos': T, 'neg': F, 'z': F}
    assert D({'pos': F}) == {'pos': F                  }
    assert D({'neg': T}) == {'pos': F, 'neg': T, 'z': F}
    assert D({'neg': F}) == {          'neg': F        }
    assert D({'z': T}) == {'pos': F, 'neg': F, 'z': T}
    assert D({'z': F}) == {                    'z': F}


def test_FactRules_deduce_multiple():
    # deduction that involves _several_ starting points
    f = FactRules(['real == pos | npos'])

    def D(facts):
        kb = FactKB(f)
        kb.deduce_all_facts(facts)
        return kb

    assert D({'real': T}) == {'real': T}
    assert D({'real': F}) == {'real': F, 'pos': F, 'npos': F}
    assert D({'pos': T}) == {'real': T, 'pos': T}
    assert D({'npos': T}) == {'real': T, 'npos': T}

    # --- key tests below ---
    assert D({'pos': F, 'npos': F}) == {'real': F, 'pos': F, 'npos': F}
    assert D({'real': T, 'pos': F}) == {'real': T, 'pos': F, 'npos': T}
    assert D({'real': T, 'npos': F}) == {'real': T, 'pos': T, 'npos': F}

    assert D({'pos': T, 'npos': F}) == {'real': T, 'pos': T, 'npos': F}
    assert D({'pos': F, 'npos': T}) == {'real': T, 'pos': F, 'npos': T}


def test_FactRules_deduce_multiple2():
    f = FactRules(['real == neg | zero | pos'])

    def D(facts):
        kb = FactKB(f)
        kb.deduce_all_facts(facts)
        return kb

    assert D({'real': T}) == {'real': T}
    assert D({'real': F}) == {'real': F, 'neg': F, 'zero': F, 'pos': F}
    assert D({'neg': T}) == {'real': T, 'neg': T}
    assert D({'zero': T}) == {'real': T, 'zero': T}
    assert D({'pos': T}) == {'real': T, 'pos': T}

    # --- key tests below ---
    assert D({'neg': F, 'zero': F, 'pos': F}) == {'real': F, 'neg': F,
             'zero': F, 'pos': F}
    assert D({'real': T, 'neg': F}) == {'real': T, 'neg': F}
    assert D({'real': T, 'zero': F}) == {'real': T, 'zero': F}
    assert D({'real': T, 'pos': F}) == {'real': T, 'pos': F}

    assert D({'real': T,           'zero': F, 'pos': F}) == {'real': T,
             'neg': T, 'zero': F, 'pos': F}
    assert D({'real': T, 'neg': F,            'pos': F}) == {'real': T,
             'neg': F, 'zero': T, 'pos': F}
    assert D({'real': T, 'neg': F, 'zero': F          }) == {'real': T,
             'neg': F, 'zero': F, 'pos': T}

    assert D({'neg': T, 'zero': F, 'pos': F}) == {'real': T, 'neg': T,
             'zero': F, 'pos': F}
    assert D({'neg': F, 'zero': T, 'pos': F}) == {'real': T, 'neg': F,
             'zero': T, 'pos': F}
    assert D({'neg': F, 'zero': F, 'pos': T}) == {'real': T, 'neg': F,
             'zero': F, 'pos': T}


def test_FactRules_deduce_base():
    # deduction that starts from base

    f = FactRules(['real  == neg | zero | pos',
                   'neg   -> real & !zero & !pos',
                   'pos   -> real & !zero & !neg'])
    base = FactKB(f)

    base.deduce_all_facts({'real': T, 'neg': F})
    assert base == {'real': T, 'neg': F}

    base.deduce_all_facts({'zero': F})
    assert base == {'real': T, 'neg': F, 'zero': F, 'pos': T}


def test_FactRules_deduce_staticext():
    # verify that static beta-extensions deduction takes place
    f = FactRules(['real  == neg | zero | pos',
                   'neg   -> real & !zero & !pos',
                   'pos   -> real & !zero & !neg',
                   'nneg  == real & !neg',
                   'npos  == real & !pos'])

    assert ('npos', True) in f.full_implications[('neg', True)]
    assert ('nneg', True) in f.full_implications[('pos', True)]
    assert ('nneg', True) in f.full_implications[('zero', True)]
    assert ('npos', True) in f.full_implications[('zero', True)]
