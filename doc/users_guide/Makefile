RST2HTML := rst2html

all:	guide.html

%.html	: %.txt
	$(RST2HTML) $< > $@

clean:
	rm -rf guide.html
