# Makefile for GDAGsim user manual

NAME=gdagsim-03
PDFFIGS=ctdag.pdf
EPSFIGS=ctdag.eps

default: gdag.pdf gdag2.ps.gz
	cp gdag.pdf gdag2.ps.gz ../$(NAME)/doc

gdag2.ps.gz: gdag2.ps
	gzip -c gdag2.ps > gdag2.ps.gz

gdag2.ps: gdag.ps
	psnup -2 -pa4 gdag.ps > gdag2.ps

gdag.ps: gdag.dvi
	dvips gdag -o gdag.ps

gdag.pdf: gdag.dvi ${PDFFIGS}
	pdflatex gdag
	pdflatex gdag

gdag.dvi: gdag.tex ${EPSFIGS}
	latex gdag
	latex gdag

%.pdf: %.eps
	epstopdf -f < $< > $@

%.eps: %.fig
	fig2dev -L ps $< > $@

clean:
	rm -f *~ *.bak *.out *.log *.dvi *.ps *.gz *.blg core gdag.pdf

full-clean:
	make clean
	rm -f *.eps *.pdf *.aux *.toc *.bbl


