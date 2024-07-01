from __future__ import print_function, division
from sympy.external import import_module


def load_ipython_formula_editor():
    from IPython.core.magic import register_line_magic
    from IPython.display import display, HTML

    html_data = """
<link rel="stylesheet" type="text/css" href="http://mathquill.com/lib/mathquill.css">
<script type="text/javascript" src="http://mathquill.com/lib/mathquill.js"></script>
<p>Expression: <span id="mathfld%i"></span></p>
<script>
function sendToIPython() {
    var kernel = IPython.notebook.kernel;
    var execcmd = 'resultVar = \"\"\"' + document.getElementById("mathfld%i").innerHTML + '\"\"\"';
    kernel.execute(execcmd);
}

function loadMathquill() {
    var script = document.createElement("script")
    script.type = "text/javascript";
    script.onload = function() {
        loadMathField();
    }
    script.src = "http://mathquill.com/lib/mathquill.js";
    document.getElementsByTagName("head")[0].appendChild(script);
}

function loadMathField() {
    var mathFieldSpan = document.getElementById('mathfld%i');


    var MQ = MathQuill.getInterface(2); // for backcompat
    var mathField = MQ.MathField(mathFieldSpan, {
      spaceBehavesLikeTab: true, // configurable
      handlers: {
    edit: function() {
        sendToIPython();
    }
      }
    });

    mathField.__controller.writeLatex("%s")
}
//loadMathquill();
loadMathField();
</script>
"""

    @register_line_magic
    def formula(expr):
        from sympy import sympify, latex
        expr = sympify(expr)
        out_html = html_data % (formula.counter, formula.counter, formula.counter, latex(expr).replace("\\", "\\\\"))
        formula.counter += 1
        display(HTML(out_html))

    formula.counter = 0
