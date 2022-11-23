PDF_FILES:=$(wildcard boundary-images/*.pdf)
SVG_FILES:=$(PDF_FILES:.pdf=.svg) 

all: boundary-variants.docx boundary-variants.pdf boundary-variants.html $(SVG_FILES)

clean:
	rm -f boundary-variants.docx boundary-variants.pdf boundary-variants.html $(SVG_FILES)

PANDOC=pandoc
PD_OPTIONS=-f markdown+implicit_figures+inline_notes+yaml_metadata_block+fenced_code_attributes --lua-filter=./pagebreak.lua --standalone
PD_EXTRA_OPTIONS=--number-sections --table-of-contents --toc-depth=1
PD_IMAGE=--default-image-extension=pdf
PD_HTML=--css ./pandoc.css --mathjax --html-q-tags
RESPATH=--resource-path=.:boundary-images

boundary-variants.docx: boundary-variants.md $(PDF_FILES)
	$(PANDOC) -o $@ -t docx $(PD_OPTIONS) $(PD_EXTRA_OPTIONS) $(PD_IMAGE) $(RESPATH) boundary-variants.md

boundary-variants.html: boundary-variants.md $(SVG_FILES)
	$(PANDOC) -o $@ -t html5 $(PD_OPTIONS) $(PD_EXTRA_OPTIONS) --default-image-extension=svg $(RESPATH) $(PD_HTML) boundary-variants.md

boundary-variants.pdf: boundary-variants.md $(PDF_FILES)
	$(PANDOC) -o $@ -t latex $(PD_OPTIONS) $(PD_EXTRA_OPTIONS) $(PD_IMAGE) $(RESPATH) boundary-variants.md

boundary-images/%.svg: boundary-images/%.pdf
	pdf2svg $< $@ 1