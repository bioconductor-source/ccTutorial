all: tex1 tex2 copycat clean

tex1: ccTutorial.tex
	pdflatex ccTutorial
	pdflatex ccTutorial

tex2: ccTutorialSupplement.tex
	pdflatex ccTutorialSupplement
	pdflatex ccTutorialSupplement

copycat: ccTutorial.tex
	cp -p ../vignettes/ccTutorial.pdf .
	cp -p ../vignettes/ccTutorialSupplement.pdf .

code:	header.tex
	R CMD Stangle ../vignettes/ccTutorial.Rnw
	R CMD Stangle ../vignettes/ccTutorialSupplement.Rnw
	cat header.tex ccTutorial.R footer.tex > codeTutorial.tex
	cat header.tex ccTutorialSupplement.R footer.tex > codeSupplement.tex
	pdflatex codeTutorial.tex
	pdflatex codeSupplement.tex
	make clean

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc
	rm -f *Tutorial*.tex code*.tex
	rm -f *.dvi *.ps ccTutorial-* ccTutorialSupplement-*
