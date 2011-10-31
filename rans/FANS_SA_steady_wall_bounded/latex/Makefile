vectorsources := $(shell find rawfigs/ ../common/rawfigs/ -name '*.dia' -o -name '*.eps' -o -name '*.ps' -o -name '*.pdf' -o -name '*.svg')
rastersources := $(shell find rawfigs/ ../common/rawfigs/ -name '*.jpg' -o -name '*.gif')
readysources := $(shell find rawfigs/ ../common/rawfigs/ -name '*.png' -o -name '*.pdf')
texsource := $(wildcard main.tex talk.tex paper.tex poster.tex proposal.tex)

vectorfigs := $(shell echo ' ' $(vectorsources) ' ' | sed -e 's> \(../common/\)*raw> >g' -e 's/\.[^. ]* /.pdf /g')
rasterfigs := $(shell echo ' ' $(rastersources) ' ' | sed -e 's> \(../common/\)*raw> >g' -e 's/\.[^. ]* /.png /g')
readyfigs := $(shell echo ' ' $(readysources) ' ' | sed -e 's> \(../common/\)*raw> >g')
figures := $(vectorfigs) $(rasterfigs) $(readyfigs)
texroot := $(patsubst %.tex, %, $(texsource))
texfinal := $(texroot).pdf

bibfiles := $(wildcard *.bib)
styfiles := $(wildcard ../common/*.sty)
clsfiles := $(wildcard ../common/*.cls)

PDFLATEX = TEXINPUTS=../common:$$TEXINPUTS pdflatex -halt-on-error
LATEX = TEXINPUTS=../common:$$TEXINPUTS latex -halt-on-error

all: $(texfinal)

figures: $(figures)

$(texfinal): *.tex $(bibfiles) $(figures) $(clsfiles) $(styfiles)
	if ls *.bib 2>&1; then $(PDFLATEX) $(texroot) && bibtex $(texroot); fi
	$(PDFLATEX) $(texroot) && $(PDFLATEX) $(texroot)

clean: cleanlatex cleanfigs

cleanlatex:
	rm -f $(patsubst %.tex, %.aux, $(wildcard *.tex)) 
	rm -f $(texroot).bbl $(texroot).blg $(texroot).log $(texroot).dvi $(texroot).nav $(texroot).out $(texroot).pdf $(texroot).snm $(texroot).toc $(texroot).vrb 

cleanfigs:
	rm -f $(figures)

figs/%.pdf: rawfigs/%.dia
	@mkdir -p $(dir $@)
	dia -t eps-builtin -e $?_roytemp.eps $? && epstopdf $?_roytemp.eps -o=$@
	@rm -f $?_roytemp.eps

figs/%.pdf: rawfigs/%.eps
	@mkdir -p $(dir $@)
	epstopdf $? -o=$@

figs/%.pdf: rawfigs/%.pdf
	@mkdir -p $(dir $@)
	@reldir=`echo $(dir $@) | sed -e 's>[^/]*/*>../>g'`; ln -sf $${reldir}$? $@

figs/%.pdf: rawfigs/%.ps
	@mkdir -p $(dir $@)
	ps2pdf $? $@

#figs/%.pdf: rawfigs/%.png
#       convert $? $@

figs/%.pdf: rawfigs/%.svg
	@mkdir -p $(dir $@)
	inkscape $? -z --export-pdf=$@

figs/%.png: rawfigs/%.jpg
	@mkdir -p $(dir $@)
	convert $? $@

figs/%.png: ../common/rawfigs/%.jpg
	@mkdir -p $(dir $@)
	convert $? $@

figs/%.png: rawfigs/%.png
	@mkdir -p $(dir $@)
	@reldir=`echo $(dir $@) | sed -e 's>[^/]*/*>../>g'`; ln -sf $${reldir}$? $@

figs/%.png: ../common/rawfigs/%.png
	@mkdir -p $(dir $@)
	@reldir=`echo $(dir $@) | sed -e 's>[^/]*/*>../>g'`; ln -sf $${reldir}$? $@

figs/%.pdf: ../common/rawfigs/%.pdf
	@mkdir -p $(dir $@)
	@reldir=`echo $(dir $@) | sed -e 's>[^/]*/*>../>g'`; ln -sf $${reldir}$? $@

figs/%.pdf: ../common/rawfigs/%.eps
	@mkdir -p $(dir $@)
	epstopdf $? -o=$@                

figs/%.pdf: ../common/rawfigs/%.ps
	@mkdir -p $(dir $@)
	ps2pdf $? $@

figs/%.pdf: ../common/rawfigs/%.svg
	@mkdir -p $(dir $@)
	inkscape $? -z --export-pdf=$@
