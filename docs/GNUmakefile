EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard lo_density/*.eps) \
             $(wildcard initial_models/*.eps) \
             $(wildcard rotation/*.eps) \
             $(wildcard radial_base/*.eps) \
             $(wildcard spherical_basestate/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = flowchart/flowchart.tex \
                  runtime_parameters/runtime_parameters.tex \
                  volume_discrepancy/volume_discrepancy.tex \
                  eos_notes/eos_notes.tex \
                  eta_notes/eta.tex \
                  pert_notes/pert.tex \
                  Godunov_notes/Godunov.tex \
                  initial_models/initial_models.tex \
                  lo_density/lo_density.tex \
                  rotation/rotation.tex \
                  radial_base/radial_base.tex \
                  spherical_basestate/basestate.tex \
                  FAQ/faq.tex \
                  thermo_notes/thermo_notes.tex

paper.dvi: paper.tex $(eps_source) $(section_sources) $(eps_source)
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null

pdf:	paper.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) paper.dvi 
	$(RM) paper.blg paper.log
	$(RM) paper.aux paper.ps paper.bbl
	$(RM) *~

.PHONY: clean
