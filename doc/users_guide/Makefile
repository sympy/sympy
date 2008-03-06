RST2HTML := rst2html

all:	guide.html

%.html	: %.txt
	$(RST2HTML) $< > $@

%.txt: src/%.txt
	cp $< $@

clean:
	rm -rf guide.html
