import os
import subprocess

here = os.path.abspath(os.path.dirname(__file__))


class Antlr4BuildError(Exception):
    pass


class LaTeXParserNotAvailable(Exception):
    pass


class Antlr4RuntimeMissing(Exception):
    pass


class Antlr4Missing(Exception):
    pass


def build_with_antlr():
    """ Attempt to call antlr4 to generate python classes from PS.g4
    """
    try:
        subprocess.check_output(["antlr4", os.path.join(here, "PS.g4")])
    except subprocess.CalledProcessError as err:
        raise Antlr4Missing("Please install antlr4")

    try:
        subprocess.check_output(["antlr4", os.path.join(here, "PS.g4")])
    except subprocess.CalledProcessError as err:
        raise Antlr4BuildError("Couldn't build LaTeX parser")


def get_parser(retry=True):
    """ Return the LaTeX parser and lexer classes
    """
    try:
        from sympy.parsing._latex.PSParser import PSParser
        from sympy.parsing._latex.PSLexer import PSLexer
        return PSParser, PSLexer
    except ModuleNotFoundError:
        print("Generated Latex parser not found, attempting to build with antlr4")
        try:
            build_with_antlr()
        except Antlr4BuildError:
            return None, None

        return get_parser(False)
