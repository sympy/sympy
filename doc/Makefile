# Makefile for Sphinx documentation
#
PYTHON       = python
RST2HTML     = rst2html

# You can set these variables from the command line.
SPHINXOPTS   =
SPHINXVER    = 0.5
#SPHINXBUILDpy  = sphinx/Sphinx-$(SPHINXVER)/sphinx-build.py
SPHINXBUILDpy  = sphinx-build
SPHINXBUILD  = PYTHONPATH=..:$(PYTHONPATH) $(SPHINXBUILDpy)
PAPER        =

SVGFILES = $(wildcard src/modules/physics/mechanics/*.svg)
PDFFILES = $(SVGFILES:%.svg=%.pdf)

ALLSPHINXOPTS = -d _build/doctrees -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) src

ALLSPHINXOPTSapi = -d _build/doctreesapi -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) api

LANGUAGES = bg cs de fr pl ru sr

.PHONY: help clean html web htmlhelp latex changes linkcheck cheatsheet

.SUFFIXES: .pdf .svg

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html        to make standalone HTML files"
	@echo "  web         to make files usable by Sphinx.web"
	@echo "  htmlhelp    to make HTML files and a HTML help project"
	@echo "  latex       to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  texinfo     to make Texinfo files"
	@echo "  info        to make Texinfo files and run them through makeinfo"
	@echo "  changes     to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck   to check all external links for integrity"
	@echo "  cheatsheet  to make the Cheatsheet"
	@echo "  htmli18n    to make the i18n versions"
	@echo "  html-errors to make the standalone HTML files, stopping on any errors or warnings"

clean:
	-rm -rf _build
	-rm -rf sphinx

$(SPHINXBUILDpy):
	rm -rf sphinx
#	mkdir sphinx
#	cd sphinx; wget http://pypi.python.org/packages/source/S/Sphinx/Sphinx-$(SPHINXVER).tar.gz;
#	cd sphinx; tar xzf Sphinx-$(SPHINXVER).tar.gz

html: $(SPHINXBUILDpy)
	mkdir -p src/.static
	mkdir -p _build/html
	mkdir -p _build/doctrees
	mkdir -p src/modules
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	cp -r src/pics _build/html/
	cp -r src/figures _build/html/
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

html-errors: SPHINXOPTS += -W
html-errors: html

htmlapi: $(SPHINXBUILDpy)
	mkdir -p api/.static
	mkdir -p api/modules
	mkdir -p _build/api _build/doctreesapi
	rm -f api/modules/sympy*.rst
	./generate_reference.py
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTSapi) _build/api
	@echo
	@echo "Build finished. The API docs pages are in _build/api."

gettext:
	mkdir -p _build/
	$(SPHINXBUILD) -b gettext src _build/gettext src/tutorial.rst
	cp _build/gettext/tutorial.pot src/

htmli18n:
	mkdir -p _build/i18n/
	mkdir -p _build/html/
	for l in $(LANGUAGES); do \
		mkdir -p src/i18n/$$l/LC_MESSAGES/ ; \
		msgfmt src/tutorial.$$l.po -o src/i18n/$$l/LC_MESSAGES/tutorial.mo ; \
		$(SPHINXBUILD) -b html -D language=$$l $(ALLSPHINXOPTS) _build/i18n/$$l src/tutorial.rst ; \
		cp _build/i18n/$$l/tutorial.html _build/html/tutorial.$$l.html ; \
	done

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

latex: $(SPHINXBUILDpy) $(PDFFILES)
	mkdir -p _build/latex _build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
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

cheatsheet:
	mkdir -p _build/cheatsheet
	pdflatex -output-directory=_build/cheatsheet cheatsheet/cheatsheet.tex
	pdflatex -output-directory=_build/cheatsheet cheatsheet/cheatsheet.tex

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
