import sys
from sympy.external import import_module

AutolevParser = AutolevLexer = AutolevListener = None
try:
    AutolevParser = import_module('sympy.parsing.autolev._antlr.autolevparser',
                                  __import__kwargs={'fromlist': ['AutolevParser']}).AutolevParser
    AutolevLexer = import_module('sympy.parsing.autolev._antlr.autolevlexer',
                                 __import__kwargs={'fromlist': ['AutolevLexer']}).AutolevLexer
    AutolevListener = import_module('sympy.parsing.autolev._antlr.autolevlistener',
                                    __import__kwargs={'fromlist': ['AutolevListener']}).AutolevListener
except Exception:
    pass


def parse_autolev(inp, output):
    antlr4 = import_module('antlr4', warn_not_installed=True)
    if not antlr4:
        raise ImportError("Autolev parsing requires the antlr4 python package,"
                          " provided by pip (antlr4-python2-runtime or"
                          " antlr4-python3-runtime) or"
                          " conda (antlr-python-runtime)")
    try:
        input_stream = antlr4.FileStream(inp)
    except Exception:
        input_stream = antlr4.InputStream(inp)

    if AutolevListener:
        from ._listener_autolev_antlr import MyListener
        lexer = AutolevLexer(input_stream)
        token_stream = antlr4.CommonTokenStream(lexer)
        parser = AutolevParser(token_stream)
        tree = parser.prog()
        my_listener = MyListener(output)
        walker = antlr4.ParseTreeWalker()
        walker.walk(my_listener, tree)
        if output == "list":
            return my_listener.output_code
