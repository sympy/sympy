"""
File to test Kamke ODEs. To run test suite, run

    python test_kamke.py

Information about failing ODEs is logged in test_report.txt.
"""


from contextlib import contextmanager
import threading
import _thread
import time

from sympy import symbols, Function, Eq, dsolve, sympify, checkodesol


class TimeOutError(Exception):
    def __init__(self, msg=''):
        self.msg = msg


@contextmanager
def time_limit(seconds, msg=''):
    timer = threading.Timer(seconds, lambda: _thread.interrupt_main())
    timer.start()
    try:
        yield
    except KeyboardInterrupt:
        raise TimeOutError("Timed out for operation {}".format(msg))
    finally:
        # Cancel timer if process finishes in time
        timer.cancel()


start = time.time()
# Read all parsed ODEs from text file
f = open("kamke_parsed.txt", "r")
eqs = f.readlines()
f.close()
# Open a file to log failing ODEs
f = open("test_report.txt", "w")

# Counts for result
solved = 0
check_slow = 0
fail = 0

# Range of ODEs to be tested
st = 0
end = len(eqs)

x = symbols('x')
y = Function('y')

for i in range(st, end):
    # Process is not terminating despite using
    # time limit. Checkodesol is too slow for
    # these cases. This should be fixed.
    if i in [13, 23, 30, 46, 47, 93, 112, 204]:
        continue

    eq = eqs[i][:-1]
    e = sympify(eq)
    # Initialize sol to a dummy value
    sol = -1
    print(i + 1, end=" ")

    try:
        # Try to find the solution to the equation
        with time_limit(10):
            sol = dsolve(e, y(x))
            if isinstance(sol, Eq):
                sol = [sol]
    except (TimeOutError, ValueError, NotImplementedError, TypeError):
        # Solution not found / timeout
        f.write(f"{i + 1} FAIL {eqs[i]}")
        fail += 1
        print("dsolve fail")

    # If a solution is found
    if sol != -1:
        try:
            # Try to verify if the solution is correct
            assert len(sol)
            with time_limit(10):
                checks = checkodesol(e, sol, y(x))
            assert any([x[0] for x in checks])
            solved += 1
            print("solved!")
        except AssertionError:
            # Wrong solution
            f.write(f"{i + 1} FAIL {eqs[i]}")
            fail += 1
            print("dsolve fail")
        except (TimeOutError, ValueError, NotImplementedError, TypeError):
            # Checkodesol unable to verify / timeout
            f.write(f"{i + 1} CHECKODESOL_FAIL {eq} {sol}/n")
            check_slow += 1
            print("checkodesol fail")

# Summary
f.write("Total time taken", time.time() - start)
f.write("No. of ODEs solved:", solved)
f.write("No. of ODEs failed:", fail)
f.write("No. of ODEs solved, but not checked:", check_slow)
f.close()

print("Total time taken", time.time() - start)
print("No. of ODEs solved:", solved)
print("No. of ODEs failed:", fail)
print("No. of ODEs solved, but not checked:", check_slow)
