from sympy.external import import_module
import os
import sys
import itertools

antlr4 = import_module("antlr4")

if not antlr4:
    disabled = True


def _test_examples(in_filename, out_filename, test_name=""):
    from sympy.parsing.autolev import parse_autolev

    dir_path = os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))
    in_file_path = os.path.join(dir_path, 'autolev', 'test-examples', in_filename)
    correct_file_path = os.path.join(dir_path, 'autolev', 'test-examples', out_filename)
    with open(in_file_path) as f:
        generated_code = parse_autolev(f, include_numeric=True)

    with open(correct_file_path) as f:
        for idx, line1 in enumerate(f):
            if line1.startswith("#"):
                break
            try:
               line2 = generated_code.split('\n')[idx]
               assert line1.rstrip() == line2.rstrip()
            except Exception:
                raise AssertionError('mismatch in ' + test_name + ' in line no: {0}'.format(idx+1))                
def test_rule_tests():
    l = ["ruletest1", "ruletest2", "ruletest3", "ruletest4", "ruletest5", "ruletest6",\
         "ruletest7", "ruletest8", "ruletest9", "ruletest10", "ruletest11", "ruletest12"]
    for i in l:
        in_filepath = i + ".al"
        out_filepath = i + ".py"
        _test_examples(in_filepath, out_filepath, i)


def test_pydy_examples():
    l = ["mass_spring_damper", "chaos_pendulum", "double_pendulum", "non_min_pendulum"]
    for i in l:
        in_filepath = os.path.join("pydy-example-repo", i + ".al")
        out_filepath = os.path.join("pydy-example-repo", i + ".py")
        _test_examples(in_filepath, out_filepath, i)


def test_autolev_tutorial():
    dir_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__)))), 'autolev', 'test-examples', 'autolev-tutorial')
    if(os.path.isdir(dir_path)):
        l = ["tutor1", "tutor2", "tutor3", "tutor4", "tutor5", "tutor6", "tutor7"]
        for i in l:
            in_filepath = os.path.join("autolev-tutorial", i + ".al")
            out_filepath = os.path.join("autolev-tutorial", i + ".py")
            _test_examples(in_filepath, out_filepath, i)


def test_dynamics_online():
    dir_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__)))), 'autolev', 'test-examples', 'dynamics-online')
    if(os.path.isdir(dir_path)):
        ch1 = ["1-4", "1-5", "1-6", "1-7", "1-8", "1-9_1", "1-9_2", "1-9_3"]
        ch2 = ["2-1", "2-2", "2-3", "2-4", "2-5", "2-6", "2-7", "2-8", "2-9", "circular"]
        ch3 = ["3-1_1", "3-1_2", "3-2_1", "3-2_2", "3-2_3", "3-2_4", "3-2_5", "3-3"]
        ch4 = ["4-1_1", "4-2_1", "4-4_1", "4-4_2", "4-5_1", "4-5_2"]
        chapters = [(ch1, "ch1"), (ch2, "ch2"), (ch3, "ch3"), (ch4, "ch4")]
        for ch, name in chapters:
            for i in ch:
                in_filepath = os.path.join("dynamics-online", name , i + ".al")
                out_filepath = os.path.join("dynamics-online", name , i + ".py")
                _test_examples(in_filepath, out_filepath, i)
