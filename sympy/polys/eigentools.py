from sympy.polys import Poly
from sympy.polys.domainmatrix import DomainMatrix
from sympy.polys.factortools import dup_factor_list
from sympy.polys.agca.extensions import FiniteExtension, ExtensionElement
from sympy.core.symbol import Dummy
from sympy.polys.rootoftools import CRootOf
from sympy.polys.polyroots import roots


def dom_eigenvects(A, l=Dummy('lambda')):
    charpoly = A.charpoly()
    rows, cols = A.shape
    domain = A.domain
    _, factors = dup_factor_list(charpoly, domain)

    eigenvects = []
    for base, exp in factors:
        if len(base) == 2:
            field = domain
            eigenval = -base[1] / base[0]

            EE_items = [
                [eigenval if i == j else field.zero for j in range(cols)]
                for i in range(rows)]
            EE = DomainMatrix(EE_items, (rows, cols), field)

            basis = (A - EE).nullspace()
            eigenvects.append((field, eigenval, exp, basis))
        else:
            minpoly = Poly.from_list(base, l, domain=domain)
            field = FiniteExtension(minpoly)
            eigenval = field(l)

            AA_items = [
                [Poly.from_list([item], l, domain=domain).rep for item in row]
                for row in A.rep]
            AA_items = [[field(item) for item in row] for row in AA_items]
            AA = DomainMatrix(AA_items, (rows, cols), field)
            EE_items = [
                [eigenval if i == j else field.zero for j in range(cols)]
                for i in range(rows)]
            EE = DomainMatrix(EE_items, (rows, cols), field)

            basis = (AA - EE).nullspace()
            eigenvects.append((field, eigenval, exp, basis))

    return eigenvects


def ddm_eigenvects_to_sympy(rational_eigenvects, Matrix, **kwargs):
    algebraic_eigenvects = []

    for field, eigenvalue, multiplicity, eigenvects in rational_eigenvects:
        if isinstance(eigenvalue, ExtensionElement):
            ring = eigenvalue.ext
            minpoly = ring.modulus
            l = minpoly.gens[0]

            eigenvects = [
                [ring.to_sympy(x) for x in vect] for vect in eigenvects]

            degree = minpoly.degree()
            minpoly = minpoly.as_expr()
            eigenvals = roots(minpoly, l, **kwargs)
            if len(eigenvals) != degree:
                eigenvals = [CRootOf(minpoly, l, idx) for idx in range(degree)]

            for eigenvalue in eigenvals:
                new_eigenvects = [
                    Matrix([x.subs(l, eigenvalue) for x in vect])
                    for vect in eigenvects]
                algebraic_eigenvects.append(
                    (eigenvalue, multiplicity, new_eigenvects))

        else:
            eigenvalue = field.to_sympy(eigenvalue)
            new_eigenvects = [
                Matrix([field.to_sympy(x) for x in vect])
                for vect in eigenvects]
            algebraic_eigenvects.append(
                (eigenvalue, multiplicity, new_eigenvects))

    return algebraic_eigenvects
