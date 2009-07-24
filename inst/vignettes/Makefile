all:	pdf

pdf:	ccTutorial.tex
	pdflatex ccTutorial
	bibtex ccTutorial
	pdflatex ccTutorial
	pdflatex ccTutorial
	pdflatex ccTutorial
	make clean

supp:	ccTutorialSupplement.tex
	pdflatex ccTutorialSupplement
	bibtex ccTutorialSupplement
	pdflatex ccTutorialSupplement
	pdflatex ccTutorialSupplement
	pdflatex ccTutorialSupplement
	make clean

supplement:	make supp

weaver:	ccTutorial.Rnw
	echo "library(weaver); Sweave(\"ccTutorial.Rnw\", driver=weaver())" | R --no-save --no-restore

nocache: echo "library(weaver); Sweave(\"ccTutorial.Rnw\", driver=weaver(), use.cache=FALSE)" | R --no-save --no-restore

weaversupp:	ccTutorialSupplement.Rnw
	echo "library(weaver); Sweave(\"ccTutorialSupplement.Rnw\", driver=weaver())" | R --no-save --no-restore

clean:
	rm -f *.out *.bbl *.aux *.log *.blg *.brf *.toc  *~ #*.eps
	rm -f ccTutorial-*.pdf ccTutorial-*.eps
	rm -f ccTutorialSupplement-*.pdf ccTutorialSupplement-*.eps
	rm -f *.nfs*
