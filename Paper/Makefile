# Define file names
ORB2_RNW_FILE=ORB2.Rnw
ORB2_TEX_FILE=ORB2.tex
ORB2_PDF_FILE=ORB2.pdf

SUPPLEMENT_RNW_FILE=Supplementary.Rnw
SUPPLEMENT_TEX_FILE=Supplementary.tex
SUPPLEMENT_PDF_FILE=Supplementary.pdf

BIB_FILE=biblio.bib

# Default target
all: clean $(ORB2_PDF_FILE) $(SUPPLEMENT_PDF_FILE)

# Rule to compile the ORB2.Rnw to ORB2.tex
$(ORB2_TEX_FILE): $(ORB2_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(ORB2_RNW_FILE)')"

# Rule to compile the ORB2.tex to ORB2.pdf
$(ORB2_PDF_FILE): $(ORB2_TEX_FILE) $(BIB_FILE)
	pdflatex $(ORB2_TEX_FILE)
	bibtex $(basename $(ORB2_TEX_FILE))
	pdflatex $(ORB2_TEX_FILE)
	pdflatex $(ORB2_TEX_FILE)

# Rule to compile the Supplement.Rnw to Supplement.tex
$(SUPPLEMENT_TEX_FILE): $(SUPPLEMENT_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(SUPPLEMENT_RNW_FILE)')"

# Rule to compile the Supplement.tex to Supplement.pdf
$(SUPPLEMENT_PDF_FILE): $(SUPPLEMENT_TEX_FILE) $(BIB_FILE)
	pdflatex $(SUPPLEMENT_TEX_FILE)
	bibtex $(basename $(SUPPLEMENT_TEX_FILE))
	pdflatex $(SUPPLEMENT_TEX_FILE)
	pdflatex $(SUPPLEMENT_TEX_FILE)

# Clean up intermediate files
clean:
	rm -f $(ORB2_TEX_FILE) $(ORB2_PDF_FILE) $(SUPPLEMENT_TEX_FILE) $(SUPPLEMENT_PDF_FILE) *.aux *.bbl *.blg *.log *.out


