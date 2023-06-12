import lark

with open('latex.lark', 'r') as f:
    latex_parser = lark.Lark(f, start='string')
    
    string = latex_parser.parse("a * b")
    print(string)
    print(string.pretty())