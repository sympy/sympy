from .utils import random_term, wilson_interval


def compute_density(presentation, variables, k=4, mode="simple"):
    return presentation.invariant_density(
        variables, k=k, mode=mode
    )


def reduce_presentation(presentation, max_depth=5):
    reduced = presentation.reduce(max_depth=max_depth)
    D0 = presentation.description_length()
    D1 = reduced.description_length()
    gain = 0 if D0 == 0 else 100 * (D0 - D1) / D0
    return reduced, gain


def estimate_invariants_mc(
    presentation,
    variables,
    k=6,
    samples=1000,
    max_depth=5,
):
    ops = [g for g in presentation.generators if callable(g)]
    hits = 0

    for _ in range(samples):
        t1 = presentation.normalize(
            random_term(variables, ops, k)
        )
        t2 = presentation.normalize(
            random_term(variables, ops, k)
        )
        if t1 != t2 and presentation.bounded_derive(
            presentation.relations[0].__class__(t1, t2),
            max_depth=max_depth,
        ):
            hits += 1

    p_hat = hits / samples
    return p_hat, wilson_interval(p_hat, samples)
