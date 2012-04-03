
from sympy import TableForm, array

def test_TableForm():
    a = array([4, 2, 3])
    s = str(TableForm(zip(a, a**3)))

    s = str(TableForm([["a", "b"], ["c", "d"], ["e", "f"]],
        headings="automatic"))
    s = str(TableForm([["a", "b"], ["c", "d"], ["e", "f"]],
            headings=("automatic", None)))
    s = str(TableForm([["a", "b"], ["c", "d"], ["e", "f"]],
            headings=(None, "automatic")))
    s = str(TableForm([[5, 7], [4, 2], [10, 3]],
            headings=[["Group A", "Group B", "Group C"], ["y1", "y2"]]))
    s = str(TableForm([[5, 7], [4, 2], [10, 3]],
            headings=[["Group A", "Group B", "Group C"], ["y1", "y2"]],
            alignment="right"))
