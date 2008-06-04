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

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	rm -f *.dvi *.ps ccTutorial-* ccTutorialSupplement-*



