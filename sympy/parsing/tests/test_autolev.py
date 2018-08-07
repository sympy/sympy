from sympy.external import import_module
import os

antlr4 = import_module("antlr4")

if not antlr4:
    disabled = True


def _test_examples(in_filename, out_filename, test_name=""):
    from sympy.parsing.autolev import parse_autolev

    dir_path = os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))
    in_file_path = os.path.join(dir_path, 'autolev', 'test_examples', in_filename)
    correct_file_path = os.path.join(dir_path, 'autolev', 'test_examples', out_filename)
    out_file_path = os.path.join(dir_path, 'autolev', 'test_examples', 'output.py')

    parse_autolev(in_file_path, out_file_path)
    #parse_autolev(in_file_path, correct_file_path)
    with open(out_file_path) as f1, open(correct_file_path) as f2:

            for idx, (lineA, lineB) in enumerate(zip(f1, f2)):
                try:
                    assert lineA.rstrip() == lineB.rstrip()
                except Exception:
                    raise AssertionError('mismatch in ' + test_name + ' in line no: {0}'.format(idx+1))


def test_rule_tests():
    number_of_tests = 12
    for i in range(1, number_of_tests+1):
        in_filepath = "ruletest" + str(i) + ".al"
        out_filepath = "ruletest" + str(i) + ".py"
        _test_examples(in_filepath, out_filepath, "test_" + str(i))

def test_mass_spring_damper():
    _test_examples("mass_spring_damper.al", "mass_spring_damper.py")

def test_chaos_pendulum():
    _test_examples("chaos_pendulum.al", "chaos_pendulum.py")

def test_double_pendulum():
    _test_examples("double_pendulum.al", "double_pendulum.py")

def test_non_min_pendulum():
    _test_examples("non_min_pendulum.al", "non_min_pendulum.py")
