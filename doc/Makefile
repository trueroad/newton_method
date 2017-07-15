all: doc

.PHONY: all doc doc-html doc-pdf clean clean-template

CSS = github-markdown.css
NUMBER_SECTIONS = --number-sections
TOC = --toc
HTML_MATH = --mathjax

CROSSREF_FILTER = --filter pandoc-crossref \
	-M "crossrefYaml=pandoc-crossref.yaml"

FILTER = $(CROSSREF_FILTER)


doc: doc-html doc-pdf
doc-html: sample-formula.ja.html sample-formula.html
doc-pdf: sample-formula.ja.pdf sample-formula.pdf


.PRECIOUS: %.bb.ja.html %.bb.ja.tex %.ja.tex \
	%.bb.html %.bb.tex %.tex \
	template.html5.org template.latex.org \
	template.ja.html5 template.ja.latex \
	template.html5 template.latex


%.bb.ja.html: %.bb.ja.md
	pandoc $< -t html5 $(HTML_MATH) -o $@

%.bb.ja.tex: %.bb.ja.md
	pandoc $< -t latex -o $@

%.bb.html: %.bb.md
	pandoc $< -t html5 $(HTML_MATH) -o $@

%.bb.tex: %.bb.md
	pandoc $< -t latex -o $@


%.ja.html: %.ja.md %.bb.ja.html template.ja.html5 $(CSS)
	pandoc $< -s $(FILTER) \
		$(NUMBER_SECTIONS) $(TOC) \
		-t html5 --template=template.ja.html5 -c $(CSS) $(HTML_MATH) \
		--include-before-body=$*.bb.ja.html \
		-o $@

%.ja.tex: %.ja.md %.bb.ja.tex template.ja.latex
	pandoc $< -s $(FILTER) \
		$(NUMBER_SECTIONS) $(TOC) \
		-t latex --template=template.ja.latex \
		-V documentclass=ltjsarticle -V colorlinks=true --listings \
		--include-before-body=$*.bb.ja.tex \
		-o $@

%.html: %.md %.bb.html template.html5 $(CSS)
	pandoc $< -s $(FILTER) \
		$(NUMBER_SECTIONS) $(TOC) \
		-t html5 --template=template.html5 -c $(CSS) $(HTML_MATH) \
		--include-before-body=$*.bb.html \
		-o $@

%.tex: %.md %.bb.tex template.latex
	pandoc $< -s $(FILTER) \
		$(NUMBER_SECTIONS) $(TOC) \
		-t latex --template=template.latex \
		-V colorlinks=true --listings \
		--include-before-body=$*.bb.tex \
		-o $@


%.pdf: %.tex
	latexmk -lualatex $^


github.css:
	wget "https://gist.github.com/andyferra/2554919/raw/10ce87fe71b23216e3075d5648b8b9e56f7758e1/github.css"

github-markdown.css:
	wget "https://github.com/sindresorhus/github-markdown-css/raw/gh-pages/github-markdown.css"


template.%.org:
	pandoc -D $* | sed -e 's/\r$$//' > $@

template.ja.%: template.ja.%.patch template.%.org
	cp template.$*.org $@
	patch -p 0 < $<

template.%: template.%.patch template.%.org
	cp template.$*.org $@
	patch -p 0 < $<


clean:
	-rm -f *~ *.tex *.pdf *.out *.toc *.aux *.log *.fdb_latexmk *.fls *.html

clean-template:
	-rm -f template.html5.org template.latex.org \
		template.ja.html5 template.ja.latex \
		template.html5 template.latex
