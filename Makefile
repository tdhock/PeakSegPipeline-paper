HOCKING-PeakSegPipeline-paper.pdf: HOCKING-PeakSegPipeline-paper.tex refs.bib
	pdflatex HOCKING-PeakSegPipeline-paper
	bibtex HOCKING-PeakSegPipeline-paper
	pdflatex HOCKING-PeakSegPipeline-paper
	pdflatex HOCKING-PeakSegPipeline-paper
figure-approx-target.png: figure-approx-target.R
	R --vanilla < $<
