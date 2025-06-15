# Define file names
#ORB2_RNW_FILE=ORB2.Rnw
#ORB2_TEX_FILE=ORB2.tex
#ORB2_PDF_FILE=ORB2.pdf

#SUPPLEMENT_RNW_FILE=Supplementary.Rnw
#SUPPLEMENT_TEX_FILE=Supplementary.tex
#SUPPLEMENT_PDF_FILE=Supplementary.pdf

#BIB_FILE=biblio.bib

# Default target
#all: clean $(ORB2_PDF_FILE) $(SUPPLEMENT_PDF_FILE)

# Rule to compile the ORB2.Rnw to ORB2.tex
#$(ORB2_TEX_FILE): $(ORB2_RNW_FILE)
#	Rscript -e "library(knitr); knitr::knit('$(ORB2_RNW_FILE)')"

# Rule to compile the ORB2.tex to ORB2.pdf
#$(ORB2_PDF_FILE): $(ORB2_TEX_FILE) $(BIB_FILE)
#	pdflatex $(ORB2_TEX_FILE)
#	bibtex $(basename $(ORB2_TEX_FILE))
#	pdflatex $(ORB2_TEX_FILE)
#	pdflatex $(ORB2_TEX_FILE)

# Rule to compile the Supplement.Rnw to Supplement.tex
#$(SUPPLEMENT_TEX_FILE): $(SUPPLEMENT_RNW_FILE)
#	Rscript -e "library(knitr); knitr::knit('$(SUPPLEMENT_RNW_FILE)')"

# Rule to compile the Supplement.tex to Supplement.pdf
#$(SUPPLEMENT_PDF_FILE): $(SUPPLEMENT_TEX_FILE)
#	pdflatex $(SUPPLEMENT_TEX_FILE)
	#bibtex $(basename $(SUPPLEMENT_TEX_FILE))
	#pdflatex $(SUPPLEMENT_TEX_FILE)
#	pdflatex $(SUPPLEMENT_TEX_FILE)

# Clean up intermediate files
#clean:
#	rm -f $(ORB2_TEX_FILE) $(ORB2_PDF_FILE) $(SUPPLEMENT_TEX_FILE) $(SUPPLEMENT_PDF_FILE) *.aux *.bbl *.blg *.log *.out

# Original files
ORB2_RNW_FILE=ORB2.Rnw
ORB2_TEX_FILE=ORB2.tex
ORB2_PDF_FILE=ORB2.pdf

SUPPLEMENT_RNW_FILE=Supplementary.Rnw
SUPPLEMENT_TEX_FILE=Supplementary.tex
SUPPLEMENT_PDF_FILE=Supplementary.pdf

# Revised files
ORB2_REVISED_RNW_FILE=ORB2_Revised.Rnw
ORB2_REVISED_TEX_FILE=ORB2_Revised.tex
ORB2_REVISED_PDF_FILE=ORB2_Revised.pdf

SUPPLEMENT_REVISED_RNW_FILE=Supplementary_Revised.Rnw
SUPPLEMENT_REVISED_TEX_FILE=Supplementary_Revised.tex
SUPPLEMENT_REVISED_PDF_FILE=Supplementary_Revised.pdf

# Diff files
ORB2_DIFF_TEX_FILE=ORB2_diff.tex
ORB2_DIFF_PDF_FILE=ORB2_diff.pdf

SUPPLEMENT_DIFF_TEX_FILE=Supplementary_diff.tex
SUPPLEMENT_DIFF_PDF_FILE=Supplementary_diff.pdf

# Response letter
LETTER_RNW_FILE=Letter.Rnw
LETTER_TEX_FILE=Letter.tex
LETTER_PDF_FILE=Letter.pdf

# Bibliography
BIB_FILE=biblio.bib

# Default target
all: clean \
     $(ORB2_PDF_FILE) $(ORB2_REVISED_PDF_FILE) $(ORB2_DIFF_PDF_FILE) \
     $(SUPPLEMENT_PDF_FILE) $(SUPPLEMENT_REVISED_PDF_FILE) $(SUPPLEMENT_DIFF_PDF_FILE) \
     $(LETTER_PDF_FILE)

# ORB2 original
$(ORB2_TEX_FILE): $(ORB2_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(ORB2_RNW_FILE)')"

$(ORB2_PDF_FILE): $(ORB2_TEX_FILE) $(BIB_FILE)
	pdflatex $(ORB2_TEX_FILE)
	bibtex $(basename $(ORB2_TEX_FILE))
	pdflatex $(ORB2_TEX_FILE)
	pdflatex $(ORB2_TEX_FILE)

# ORB2 revised
$(ORB2_REVISED_TEX_FILE): $(ORB2_REVISED_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(ORB2_REVISED_RNW_FILE)')"

$(ORB2_REVISED_PDF_FILE): $(ORB2_REVISED_TEX_FILE) $(BIB_FILE)
	pdflatex $(ORB2_REVISED_TEX_FILE)
	bibtex $(basename $(ORB2_REVISED_TEX_FILE))
	pdflatex $(ORB2_REVISED_TEX_FILE)
	pdflatex $(ORB2_REVISED_TEX_FILE)

# ORB2 diff
$(ORB2_DIFF_TEX_FILE): $(ORB2_TEX_FILE) $(ORB2_REVISED_TEX_FILE)
	latexdiff $(ORB2_TEX_FILE) $(ORB2_REVISED_TEX_FILE) > $(ORB2_DIFF_TEX_FILE)

$(ORB2_DIFF_PDF_FILE): $(ORB2_DIFF_TEX_FILE)
	pdflatex $(ORB2_DIFF_TEX_FILE)
	pdflatex $(ORB2_DIFF_TEX_FILE)

# Supplement original
$(SUPPLEMENT_TEX_FILE): $(SUPPLEMENT_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(SUPPLEMENT_RNW_FILE)')"

$(SUPPLEMENT_PDF_FILE): $(SUPPLEMENT_TEX_FILE)
	pdflatex $(SUPPLEMENT_TEX_FILE)
	pdflatex $(SUPPLEMENT_TEX_FILE)

# Supplement revised
$(SUPPLEMENT_REVISED_TEX_FILE): $(SUPPLEMENT_REVISED_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(SUPPLEMENT_REVISED_RNW_FILE)')"

$(SUPPLEMENT_REVISED_PDF_FILE): $(SUPPLEMENT_REVISED_TEX_FILE)
	pdflatex $(SUPPLEMENT_REVISED_TEX_FILE)
	pdflatex $(SUPPLEMENT_REVISED_TEX_FILE)

# Supplement diff
$(SUPPLEMENT_DIFF_TEX_FILE): $(SUPPLEMENT_TEX_FILE) $(SUPPLEMENT_REVISED_TEX_FILE)
	latexdiff $(SUPPLEMENT_TEX_FILE) $(SUPPLEMENT_REVISED_TEX_FILE) > $(SUPPLEMENT_DIFF_TEX_FILE)

$(SUPPLEMENT_DIFF_PDF_FILE): $(SUPPLEMENT_DIFF_TEX_FILE)
	pdflatex $(SUPPLEMENT_DIFF_TEX_FILE)
	pdflatex $(SUPPLEMENT_DIFF_TEX_FILE)

# Response letter
$(LETTER_TEX_FILE): $(LETTER_RNW_FILE)
	Rscript -e "library(knitr); knitr::knit('$(LETTER_RNW_FILE)')"

$(LETTER_PDF_FILE): $(LETTER_TEX_FILE)
	pdflatex $(LETTER_TEX_FILE)
	pdflatex $(LETTER_TEX_FILE)

# Clean up intermediate and output files
clean:
	rm -f \
		$(ORB2_TEX_FILE) $(ORB2_PDF_FILE) \
		$(ORB2_REVISED_TEX_FILE) $(ORB2_REVISED_PDF_FILE) \
		$(ORB2_DIFF_TEX_FILE) $(ORB2_DIFF_PDF_FILE) \
		$(SUPPLEMENT_TEX_FILE) $(SUPPLEMENT_PDF_FILE) \
		$(SUPPLEMENT_REVISED_TEX_FILE) $(SUPPLEMENT_REVISED_PDF_FILE) \
		$(SUPPLEMENT_DIFF_TEX_FILE) $(SUPPLEMENT_DIFF_PDF_FILE) \
		$(LETTER_TEX_FILE) $(LETTER_PDF_FILE) \
		*.aux *.bbl *.blg *.log *.out *.toc


