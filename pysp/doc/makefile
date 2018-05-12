HTML=doc.html

all: $(HTML)

clean:
	rm -f $(HTML)

%.html: %.txt
	asciidoc $(<F) || rm $(HTML)

.PHONY: all clean 
