from __future__ import annotations
from sympy.combinatorics.rewritingsystem_fsm import (
    epsilon_closure, instructions,
    subset_construction, subset_construction_epsilon)


def test_doctest_copies() -> None:
    transition: dict[int, dict[str, set[int]]]  = {
        0: {'e': {1, 7}},
        1: {'e': {2, 4}},
        2: {'a': {3}},
        3: {'e': {6}},
        4: {'b': {5}},
        5: {'e': {6}},
        6: {'e': {1, 7}},
        7: {'a': {8}},
        8: {'b': {9}},
        9: {'b': {10}},
    }
    epsilon = 'e'

    assert epsilon_closure(transition, epsilon, {0}) == {0, 1, 2, 4, 7}
    assert epsilon_closure(transition, epsilon, {5}) == {1, 2, 4, 5, 6, 7}
    assert epsilon_closure(transition, epsilon, {3, 8}) == {1, 2, 3, 4, 6, 7, 8}
    assert instructions(transition, {0, 1, 2, 4, 7}) == {
        'a': {8, 3}, 'b': {5}, 'e': {1, 2, 4, 7}}

    assert subset_construction(transition, {0}) == {
        frozenset({0}): {'e': frozenset({1, 7})},
        frozenset({1, 7}): {'e': frozenset({2, 4}), 'a': frozenset({8})},
        frozenset({8}): {'b': frozenset({9})},
        frozenset({9}): {'b': frozenset({10})},
        frozenset({2, 4}): {'a': frozenset({3}), 'b': frozenset({5})},
        frozenset({5}): {'e': frozenset({6})},
        frozenset({6}): {'e': frozenset({1, 7})},
        frozenset({3}): {'e': frozenset({6})}
    }
    assert subset_construction_epsilon(
        transition, 'e', epsilon_closure(transition, 'e', {0})
    ) == {
        frozenset({0, 1, 2, 4, 7}): {
            'a': frozenset({1, 2, 3, 4, 6, 7, 8}),
            'b': frozenset({1, 2, 4, 5, 6, 7})},
        frozenset({1, 2, 4, 5, 6, 7}): {
            'a': frozenset({1, 2, 3, 4, 6, 7, 8}),
            'b': frozenset({1, 2, 4, 5, 6, 7})},
        frozenset({1, 2, 3, 4, 6, 7, 8}): {
            'a': frozenset({1, 2, 3, 4, 6, 7, 8}),
            'b': frozenset({1, 2, 4, 5, 6, 7, 9})},
        frozenset({1, 2, 4, 5, 6, 7, 9}): {
            'a': frozenset({1, 2, 3, 4, 6, 7, 8}),
            'b': frozenset({1, 2, 4, 5, 6, 7, 10})},
        frozenset({1, 2, 4, 5, 6, 7, 10}): {
            'a': frozenset({1, 2, 3, 4, 6, 7, 8}),
            'b': frozenset({1, 2, 4, 5, 6, 7})
        }
    }
