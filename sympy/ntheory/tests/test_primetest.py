from sympy import Sieve
from sympy.core.compatibility import range

from sympy.ntheory import isprime
from sympy.ntheory.primetest import _mr_safe_helper


def test_isprime():
    s = Sieve()
    s.extend(100000)
    ps = set(s.primerange(2, 100001))
    for n in range(100001):
        # if (n in ps) != isprime(n): print n
        assert (n in ps) == isprime(n)
    assert isprime(179424673)
    # Some Mersenne primes
    assert isprime(2**61 - 1)
    assert isprime(2**89 - 1)
    assert isprime(2**607 - 1)
    assert not isprime(2**601 - 1)
    #Arnault's number
    assert isprime(int('''
803837457453639491257079614341942108138837688287558145837488917522297\
427376533365218650233616396004545791504202360320876656996676098728404\
396540823292873879185086916685732826776177102938969773947016708230428\
687109997439976544144845341155872450633409279022275296229414984230688\
1685404326457534018329786111298960644845216191652872597534901'''))
    # pseudoprime that passes the base set [2, 3, 7, 61, 24251]
    assert not isprime(9188353522314541)

    assert _mr_safe_helper(
        "if n < 170584961: return mr(n, [350, 3958281543])") == \
        ' # [350, 3958281543] stot = 1 clear [2, 3, 5, 7, 29, 67, 679067]'
    assert _mr_safe_helper(
        "if n < 3474749660383: return mr(n, [2, 3, 5, 7, 11, 13])") == \
        ' # [2, 3, 5, 7, 11, 13] stot = 7 clear == bases'
    assert isprime(5.0) is False
