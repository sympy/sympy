import glob

# Mapping of broken URLs to replacements
replacements = {
    "https://encyclopediaofmath.org/wiki/Airy_functions": "https://mathworld.wolfram.com/AiryFunctions.html",
    "https://github.com/sympy/sympy/blob/master/sympy/parsing/autolev/test-examples/pydy-example-repo/mass_spring_damper.py#L10":
        "https://github.com/sympy/sympy/blob/master/sympy/parsing/autolev/test-examples/pydy-example-repo/mass_spring_damper.py",
    "https://github.com/sympy/sympy/issues/21177#issuecomment-812816346": "https://github.com/sympy/sympy/issues/21177",
    "https://groupprops.subwiki.org/wiki/Structure_theorem_for_finitely_generated_abelian_groups":
        "https://en.wikipedia.org/wiki/Structure_theorem_for_finitely_generated_abelian_groups",
    "https://groupprops.subwiki.org/wiki/Derived_subgroup": "https://en.wikipedia.org/wiki/Derived_subgroup",
    "https://math24.net/bessel-differential-equation.html": "https://en.wikipedia.org/wiki/Bessel_differential_equation",
    "https://mae.ufl.edu/~fregly/PDFs/autolev_tutorial.pdf":
        "https://web.archive.org/web/20180731093609/http://docs.sympy.org/0.7.6/modules/mpmath/calculus/polynomials.html",
    "https://mersenneforum.org/showpost.php?p=110896": "https://www.mersenneforum.org/",
    "https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470316887.app1": "https://onlinelibrary.wiley.com/doi/10.1002/9780470316887.app1",
    "https://lpsa.swarthmore.edu/Representations/SysRepSS.html#:~:text=Example%3A%20Direct%20Derivation%20of%20State%20Space%20Model%20(Electrical),-Derive%20a%20state&text=The%20input%20is%20ia%20and%20the%20output%20is%20e2.&text=space%20representation%20becomes-,This%20technique%20does%20not%20always%20easily%20yield%20a%20set%20of,Transfer%20functions%20are%20discussed%20elsewhere.":
        "https://lpsa.swarthmore.edu/Representations/SysRepSS.html",
    "https://www.math.usm.edu/perry/Research/Thesis_DRL.pdf":
        "https://web.archive.org/web/20221029115428/https://web.cs.elte.hu/~lovasz/scans/lll.pdf",
    "https://pubs.aip.org/aapt/ajp/article-abstract/20/1/1/1034555/Units-and-Dimensions-in-Physics":
        "https://www.nist.gov/pml/owm/metric-si",
    "https://www.circuitbread.com/tutorials/routh-hurwitz-criterion-part-2-3-3":
        "https://en.wikipedia.org/wiki/Routh–Hurwitz_stability_criterion",
    "https://xy-pic.sourceforge.net/": "https://ctan.org/pkg/xypic",
    "https://people.brandeis.edu/~igusa/Math56aS08/Math56a_S08_notes015.pdf": "",
    "https://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks": "https://docs.python.org/3/tutorial/index.html",
    "https://mathworld.wolfram.com/DirichletDistribution.html": "https://en.wikipedia.org/wiki/Dirichlet_distribution",
}

# File patterns to search
file_patterns = ["**/*.rst", "**/*.md"]

for pattern in file_patterns:
    for filename in glob.glob(pattern, recursive=True):
        with open(filename, "r", encoding="utf-8") as f:
            content = f.read()
        new_content = content
        for old, new in replacements.items():
            new_content = new_content.replace(old, new)
        if new_content != content:
            print(f"Updated links in {filename}")
            with open(filename, "w", encoding="utf-8") as f:
                f.write(new_content)