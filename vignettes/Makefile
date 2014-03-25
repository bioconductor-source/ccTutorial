all: buildclean

%.pdf: %.tex
	cp ../inst/vignettes/$@ ./

buildclean: ccTutorial.pdf ccTutorialSupplement.pdf
	rm -f *.tex
