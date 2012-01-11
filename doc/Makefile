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

ALLSPHINXOPTS = -d _build/doctrees -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) src

ALLSPHINXOPTSapi = -d _build/doctreesapi -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) api

.PHONY: help clean html web htmlhelp latex changes linkcheck cheatsheet

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  web        to make files usable by Sphinx.web"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes    to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  cheatsheet to make the Cheatsheet"

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

htmlapi: $(SPHINXBUILDpy)
	mkdir -p api/.static
	mkdir -p api/modules
	mkdir -p _build/api _build/doctreesapi
	rm -f api/modules/sympy*.txt
	./generate_reference.py
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTSapi) _build/api
	@echo
	@echo "Build finished. The API docs pages are in _build/api."

gettext:
	mkdir -p _build/
	$(SPHINXBUILD) -b gettext src _build/gettext src/tutorial.txt
	cp _build/gettext/tutorial.pot src/

htmli18n:
	mkdir -p _build/i18n/
	mkdir -p _build/html/
	mkdir -p src/i18n/cs/LC_MESSAGES/
	mkdir -p src/i18n/ru/LC_MESSAGES/
	mkdir -p src/i18n/bg/LC_MESSAGES/
	msgfmt src/tutorial.cs.po -o src/i18n/cs/LC_MESSAGES/tutorial.mo
	msgfmt src/tutorial.ru.po -o src/i18n/ru/LC_MESSAGES/tutorial.mo
	msgfmt src/tutorial.bg.po -o src/i18n/bg/LC_MESSAGES/tutorial.mo
	$(SPHINXBUILD) -b html -D language=cs $(ALLSPHINXOPTS) _build/i18n/cs src/tutorial.txt
	$(SPHINXBUILD) -b html -D language=ru $(ALLSPHINXOPTS) _build/i18n/ru src/tutorial.txt
	$(SPHINXBUILD) -b html -D language=bg $(ALLSPHINXOPTS) _build/i18n/bg src/tutorial.txt
	cp _build/i18n/cs/tutorial.html _build/html/tutorial.cs.html
	cp _build/i18n/ru/tutorial.html _build/html/tutorial.ru.html
	cp _build/i18n/bg/tutorial.html _build/html/tutorial.bg.html

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

latex: $(SPHINXBUILDpy)
	mkdir -p _build/latex _build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

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
