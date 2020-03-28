def gauss(A_mat, modulus):
    import copy
    A =  copy.deepcopy(A_mat)
    from sympy import mod_inverse
    n = len(A)

    for i in range(n):
        for j in range(len(A[i])):
            A[i][j] =  A[i][j]  % modulus
    n_c = len(A[0])

    for i in range(0, n_c - 1):
        maxElement = A[i][i]
        max_r = i
        for k in range(i + 1, n):
            if A[k][i] > maxElement:
                maxElement = A[k][i]
                max_r = k

        for k in range(i,  n_c):
            tmp = A[max_r][k]
            A[max_r][k] = A[i][k]
            A[i][k] = tmp
        for k in range(i + 1, n):
            m = -A[k][i] * mod_inverse(A[i][i], modulus) % modulus
            for j in range(i, n_c):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j]  =  (A[k][j] + m * A[i][j])% modulus

    A = A[:n_c - 1]
    n = len(A)
    sol = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        sol[i] = A[i][n]  * mod_inverse(A[i][i], modulus) % modulus
        for k in range(i - 1, -1, -1):
            A[k][n] = (A[k][n] - A[k][i] * sol[i]) % modulus
    return sol

def is_B_smooth(n, B):
    from sympy import sieve
    b = []
    for i in  sieve.primerange(1, B + 1):
        count = 0
        while (n % i == 0):
            n //= i
            count += 1
        b.append(count)
    return n == 1, b

def index_calculus(g, base, p, B=5 , c=10):
    """
    Index calculus algorithm for discrete
    log problem.

    References
    ==========

    .. [1] Jeffrey Hoffstein
        "An introduction to mathematical cryptography" 162-165
    """
    from sympy import sieve, mod_inverse, primefactors, isprime
    import random
    from sympy.ntheory.modular import solve_congruence
    assert isprime(p )
    Bound = [i for i in sieve.primerange(1,B + 1)]
    total_independent = 0
    final_A = []
    while total_independent < len(Bound) + c:
        exp = random.randint(1, p)
        is_B, ls = is_B_smooth(pow(g, exp, p), B)
        if not is_B:
            continue
        ls.append(exp)
        final_A.append(ls)
        total_independent += 1

    pp = primefactors(p - 1)
    sols = []
    pri =[]
    for i in pp:
        sols.append(gauss(final_A, i))
        pri.append(i)
    final_sols= []
    for i in range(len(sols[0])):
        s = []
        for j in range(len(sols)):
            s.append(sols[j][i])
        final_sols.append(solve_congruence(*zip(s, pri)))
    while(1):
        a = random.randint(1,p)
        a1 = base * pow( mod_inverse(g,p) , a , p) % p
        is_B, ls = is_B_smooth(a1, B)
        if is_B:
            break

    solution = a
    for idx, j in enumerate(ls):
        solution +=   ls[idx] * final_sols[idx][0]
    return solution % (p - 1)
