# Define file names
RNW_FILE=ORB2.Rnw
TEX_FILE=ORB2.tex
PDF_FILE=ORB2.pdf
BIB_FILE=biblio.bib

# Default target
all: $(PDF_FILE)

# Rule to compile the .Rnw to .tex
$(TEX_FILE): $(RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(RNW_FILE)')"

# Rule to compile the .tex to .pdf
$(PDF_FILE): $(TEX_FILE) $(BIB_FILE)
	pdflatex $(TEX_FILE)
	bibtex ORB2
	pdflatex $(TEX_FILE)
	pdflatex $(TEX_FILE)

# Clean up intermediate files
clean:
	rm -f $(TEX_FILE) $(PDF_FILE) ORB2.aux ORB2.bbl ORB2.blg ORB2.log ORB2.out

