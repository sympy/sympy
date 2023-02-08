from collections import Counter, defaultdict

from sympy.tensor.array.expressions import ArrayTensorProduct, ArrayContraction, ArrayDiagonal, PermuteDims


def einsum_to_sympy_array(path, *args):
    path = path.replace(" ", "").replace("\t", "").strip()
    indices1, indices2 = path.split("->")
    indices1_by_arg = indices1.split(",")
    assert len(indices1_by_arg) == len(args)
    base = ArrayTensorProduct(*args)
    counted_start = Counter([j for i in indices1_by_arg for j in i])
    counted_end = Counter([i for i in indices2])
    indices_contraction = {i for i, count in counted_start.items() if count > 1 and i not in counted_end}
    indices_diagonalization = [i for i, count in counted_start.items() if count > 1 and i in counted_end]
    start_indices_pos = defaultdict(list)
    indices1_flat = [j for i in indices1_by_arg for j in i]
    order_after_contraction = [i for i in indices1_flat if i not in indices_contraction]
    order_after_diagonalization = [i for i in order_after_contraction if i not in indices_diagonalization]
    order_after_diagonalization_part2 = [i for i in order_after_contraction if i in indices_diagonalization]
    [order_after_diagonalization.append(i) for i in order_after_diagonalization_part2 if i not in order_after_diagonalization]
    for i, ind in enumerate(indices1_flat):
        start_indices_pos[ind].append(i)
    contraction_indices = []
    for i in indices_contraction:
        contraction_indices.append(start_indices_pos[i])
    diagonalization_indices = []
    for i in indices_diagonalization:
        diagonalization_indices.append(start_indices_pos[i])
    array_contraction = ArrayContraction(base, *contraction_indices)
    diagonalization_indices = array_contraction._push_indices_up(
        array_contraction.contraction_indices, diagonalization_indices)
    array_diagonal = ArrayDiagonal(array_contraction, *diagonalization_indices)
    permu = PermuteDims(array_diagonal, index_order_old=order_after_diagonalization, index_order_new=indices2)
    return permu
