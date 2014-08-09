# Makefile for Sphinx documentation
#
PYTHON       = python
RST2HTML     = rst2html

# You can set these variables from the command line.
SPHINXOPTS   =
SPHINXVER    = 0.5
SPHINXBUILD  = PYTHONPATH=..:$(PYTHONPATH) sphinx-build
PAPER        =

SVGFILES = $(wildcard src/modules/physics/vector/*.svg) $(wildcard src/modules/physics/mechanics/examples/*.svg)
PDFFILES = $(SVGFILES:%.svg=%.pdf)

ALLSPHINXOPTS = -d _build/doctrees $(SPHINXOPTS) src
ALLSPHINXOPTSapi = -d _build/doctrees-api $(SPHINXOPTS) api
ALLSPHINXOPTSlatex = -d _build/doctrees-latex -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) src

.PHONY: changes cheatsheet clean help html html-errors \
        htmlapi htmlhelp info latex linkcheck livehtml \
        texinfo web logo man

.SUFFIXES: .pdf .svg

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  changes     to make an overview over all changed/added/deprecated items"
	@echo "  cheatsheet  to make the Cheatsheet"
	@echo "  clean       to remove generated files"
	@echo "  html        to make standalone HTML files"
	@echo "  html-errors to make the standalone HTML files, stopping on any errors or warnings"
	@echo "  htmlapi     to make HTML API docs"
	@echo "  htmlhelp    to make HTML files and a HTML help project"
	@echo "  info        to make Texinfo files and run them through makeinfo"
	@echo "  latex       to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  linkcheck   to check all external links for integrity"
	@echo "  livehtml    to use livereload to view the built html"
	@echo "  texinfo     to make Texinfo files"
	@echo "  web         to make files usable by Sphinx.web"
	@echo "  logo        generate logos from src/logo/sympy.svg"
	@echo "  man         build manpage"

clean:
	-rm -rf _build
	-rm -rf sphinx
	-rm -f $(PDFFILES)

html: _build/logo/sympy-notailtext-favicon.ico
	mkdir -p src/.static
	mkdir -p _build/html
	mkdir -p _build/doctrees
	mkdir -p src/modules
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	cp -r src/pics _build/html/
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

html-errors: SPHINXOPTS += -W
html-errors: html

htmlapi:
	mkdir -p api/.static
	mkdir -p api/modules
	mkdir -p _build/api _build/doctreesapi
	rm -f api/modules/sympy*.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTSapi) _build/api
	@echo
	@echo "Build finished. The API docs pages are in _build/api."

web:
	mkdir -p _build/web _build/doctrees
	$(SPHINXBUILD) -b web $(ALLSPHINXOPTS) _build/web
	@echo
	@echo "Build finished; now you can run"
	@echo "  python -m sphinx.web _build/web"
	@echo "to start the server."

htmlhelp:
	mkdir -p _build/htmlhelp _build/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

latex: $(PDFFILES)
	mkdir -p _build/latex _build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTSlatex) _build/latex
	sed -i "s/pdflatex/xelatex/g" _build/latex/Makefile
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all' in that directory to run these through xelatex."

.svg.pdf:
	inkscape --file=$< --export-area-drawing --without-gui --export-pdf=$@

changes:
	mkdir -p _build/changes _build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	mkdir -p _build/linkcheck _build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

livehtml:
	livereload _build/html

cheatsheet: _build/cheatsheet/cheatsheet.pdf _build/cheatsheet/combinatoric_cheatsheet.pdf

_build/cheatsheet/cheatsheet.pdf: cheatsheet/cheatsheet.tex
	mkdir -p _build/cheatsheet
	pdflatex -output-directory=_build/cheatsheet cheatsheet/cheatsheet.tex
	pdflatex -output-directory=_build/cheatsheet cheatsheet/cheatsheet.tex

_build/cheatsheet/cheatsheet.pdf: cheatsheet/combinatoric_cheatsheet.tex
	mkdir -p _build/cheatsheet
	pdflatex -output-directory=_build/cheatsheet cheatsheet/combinatoric_cheatsheet.tex
	pdflatex -output-directory=_build/cheatsheet cheatsheet/combinatoric_cheatsheet.tex

texinfo:
	mkdir -p _build/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) _build/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in _build/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	mkdir -p _build/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) _build/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C _build/texinfo info
	@echo "makeinfo finished; the Info files are in _build/texinfo."

man: man/isympy.xml
	docbook2x-man --to-stdout $< > man/isympy.1

_build/logo/sympy-notailtext-favicon.ico: logo

logo: src/logo/sympy.svg
	rm -rf _build/logo
	mkdir -p _build/logo
	$(PYTHON) ./generate_logos.py -d
	@echo
	@echo "Logo generated."
