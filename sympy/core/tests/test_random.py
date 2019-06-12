import random
from sympy import Symbol, symbols, sin, acos
from sympy.abc import x as xx


# Values should be equal:
random.seed(42)
before = random.random()

random.seed(42)
Symbol('z').is_finite
after = random.random()
assert before == after

# --------------------

results = set()
for _ in range(10):
    random.seed(42)
    (xx**2).as_numer_denom()
    results.add(random.random())
# This should be 1:
assert len(results) == 1

# --------------------

# Values should be equal:
random.seed(28)
m0, m1 = symbols('m_0 m_1', real=True)
result = acos(-m0/m1)
before = random.uniform(0,1)

random.seed(28)
m0, m1 = symbols('m_0 m_1', real=True)
result = acos(-m0/m1)
after = random.uniform(0,1)

assert before == after

# --------------------

# Values should be equal:
random.seed(10)
x = Symbol('x')
y = 0
for i in range(4): y += sin(random.uniform(-10,10) * x)
before = y

random.seed(10)
x = Symbol('x')
z = 0
for i in range(4): z += sin(random.uniform(-10,10) * x)
after = z

assert before == after
