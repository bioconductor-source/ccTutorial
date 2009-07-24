all:	pdf supp

pdf:	ccTutorial.tex
	cp ../vignettes/ccTutorial.pdf .
	R CMD Stangle ../vignettes/ccTutorial.Rnw
	cp ../vignettes/ccTutorial.R .

supp:	ccTutorialSupplement.tex
	cp ../vignettes/ccTutorialSupplement.pdf .
	R CMD Stangle ../vignettes/ccTutorialSupplement.Rnw
	cp ../vignettes/ccTutorialSupplement.R .
