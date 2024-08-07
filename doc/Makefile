# Makefile for Sphinx documentation
#
PYTHON       = python
RST2HTML     = rst2html

# You can set these variables from the command line.
SPHINXOPTS   =
SPHINXVER    = 0.5
SPHINXBUILD  = sphinx-build
# This one can also be set in the environment
PAPER        ?= letter
BUILDDIR     = _build
SOURCEDIR    = src
LIVEHOST     = localhost
LIVEPORT     = 0

SVGFILES = $(wildcard src/modules/physics/vector/*.svg) $(wildcard src/modules/physics/mechanics/examples/*.svg) $(wildcard src/modules/vector/*.svg)
PDFFILES = $(SVGFILES:%.svg=%.pdf)

ALLSPHINXOPTS = -d $(BUILDDIR)/doctrees $(SPHINXOPTS)
WARNINGSFILE = $(BUILDDIR)/.sphinx-warnings
SPHINXSAVEWARNINGS = -w $(WARNINGSFILE)
ALLSPHINXOPTSapi = -d $(BUILDDIR)/doctrees-api $(SPHINXOPTS) api
ALLSPHINXOPTSlatex = -d $(BUILDDIR)/doctrees-latex -D latex_elements.papersize=$(PAPER)paper \
                $(SPHINXOPTS) src

.PHONY: changes cheatsheet clean help html _html htmlapi htmlhelp info latex \
        pdflatex latexpdf pdf linkcheck livehtml texinfo web logo man printwarnings

.SUFFIXES: .pdf .svg

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  changes     to make an overview over all changed/added/deprecated items"
	@echo "  cheatsheet  to make the Cheatsheet"
	@echo "  clean       to remove generated files"
	@echo "  html        to make standalone HTML files"
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
	-rm -rf $(BUILDDIR)
	-rm -rf sphinx
	-rm -f $(PDFFILES)

html:
	$(MAKE) -s _html || $(MAKE) -s printwarnings

# Force colors on CI. Set NO_COLOR=1 on the command line to disable this.
_html: export FORCE_COLOR = 1
_html: $(BUILDDIR)/html/pics/*.png
_html: SPHINXOPTS += -W --keep-going
_html: $(BUILDDIR)/logo/sympy-notailtext-favicon.ico
	mkdir -p $(SOURCEDIR)/.static
	mkdir -p $(BUILDDIR)/html
	mkdir -p $(BUILDDIR)/doctrees
	mkdir -p $(SOURCEDIR)/modules
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(SPHINXSAVEWARNINGS) $(SOURCEDIR) $(BUILDDIR)/html
	@echo
	@echo "Docs build finished. The HTML pages are in $(BUILDDIR)/html."

# Print the warnings again at the end of the build, to make them easier to find
printwarnings:
	@printf "\033[0;31mThe docs build failed with the following errors:\033[0m\n"
	@cat $(WARNINGSFILE)
	exit 1

$(BUILDDIR)/html/pics/*.png: $(SOURCEDIR)/pics/*.png
	mkdir -p $(BUILDDIR)/html
	cp -r $(SOURCEDIR)/pics $(BUILDDIR)/html/

htmlapi:
	mkdir -p api/.static
	mkdir -p api/modules
	mkdir -p $(BUILDDIR)/api $(BUILDDIR)/doctreesapi
	rm -f api/modules/sympy*.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTSapi) $(BUILDDIR)/api
	@echo
	@echo "Build finished. The API docs pages are in $(BUILDDIR)/api."

web:
	mkdir -p $(BUILDDIR)/web $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b web $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/web
	@echo
	@echo "Build finished; now you can run"
	@echo "  python -m sphinx.web $(BUILDDIR)/web"
	@echo "to start the server."

htmlhelp:
	mkdir -p $(BUILDDIR)/htmlhelp $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) %(SOURCEDIR) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

latex: logo
	mkdir -p $(BUILDDIR)/latex $(BUILDDIR)/doctrees-latex
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTSlatex) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make all' in that directory to run these through xelatex."

pdf pdflatex latexpdf: latex
	@$(MAKE) -C $(BUILDDIR)/latex all-pdf \
	LATEXOPTS="-halt-on-error -interaction=batchmode" \
	LATEXMKOPTS="-xelatex" \
	|| echo "LaTeX error"'!'" Check .log in $(BUILDDIR)/latex for lineno in .tex file."

changes:
	mkdir -p $(BUILDDIR)/changes $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	mkdir -p $(BUILDDIR)/linkcheck $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

livehtml: _html
	sphinx-autobuild --open-browser --watch .. --host $(LIVEHOST) --port $(LIVEPORT) -b html $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html

cheatsheet: $(BUILDDIR)/cheatsheet/cheatsheet.pdf $(BUILDDIR)/cheatsheet/combinatoric_cheatsheet.pdf

$(BUILDDIR)/cheatsheet/cheatsheet.pdf: cheatsheet/cheatsheet.tex
	mkdir -p $(BUILDDIR)/cheatsheet
	pdflatex -output-directory=$(BUILDDIR)/cheatsheet cheatsheet/cheatsheet.tex
	pdflatex -output-directory=$(BUILDDIR)/cheatsheet cheatsheet/cheatsheet.tex

$(BUILDDIR)/cheatsheet/combinatoric_cheatsheet.pdf: cheatsheet/combinatoric_cheatsheet.tex
	mkdir -p $(BUILDDIR)/cheatsheet
	pdflatex -output-directory=$(BUILDDIR)/cheatsheet cheatsheet/combinatoric_cheatsheet.tex
	pdflatex -output-directory=$(BUILDDIR)/cheatsheet cheatsheet/combinatoric_cheatsheet.tex

texinfo:
	mkdir -p $(BUILDDIR)/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	mkdir -p $(BUILDDIR)/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

man: man/isympy.xml
	docbook2x-man --to-stdout $< > man/isympy.1

logo: $(BUILDDIR)/logo/*.png

$(BUILDDIR)/logo/sympy-notailtext-favicon.ico: logo

$(BUILDDIR)/logo/*.png: $(SOURCEDIR)/logo/sympy.svg
	mkdir -p $(BUILDDIR)/logo
	$(PYTHON) ./generate_logos.py -d
	@echo
	@echo "Logo generated."

autobuild:
	sphinx-autobuild $(ALLSPHINXOPTS) $(SPHINXSAVEWARNINGS) $(SOURCEDIR) $(BUILDDIR)/html
