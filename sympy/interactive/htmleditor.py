from IPython.core.magic import register_line_magic
from IPython.display import display, HTML

html_data = """
<link rel="stylesheet" type="text/css" href="http://mathquill.com/lib/mathquill.css">
<script type="text/javascript" src="http://mathquill.com/lib/mathquill.js"></script>
<p>Expression: <span id="mathfld%i"></span>
<script>
    var mathFieldSpan = document.getElementById('mathfld%i');

    var MQ = MathQuill.getInterface(2); // for backcompat
    var mathField = MQ.MathField(mathFieldSpan, {
      spaceBehavesLikeTab: true, // configurable
      handlers: {
        edit: function() { // useful event handlers
          //latexSpan.textContent = mathField.latex(); // simple API
        }
      }
    });
</script>
<script>
mathField.__controller.writeLatex("%s")
</script>
"""

formula_counter = 0


@register_line_magic
def formula(expr):
    from sympy import sympify, latex
    global formula_counter
    expr = sympify(expr)
    out_html = html_data % (formula_counter, formula_counter, latex(expr).replace("\\", "\\\\"))
    formula_counter += 1
    display(HTML(out_html))


del formula
