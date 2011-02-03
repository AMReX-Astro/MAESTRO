EPStoPDF = epstopdf

ALL: MaestroUsersGuide.ps MaestroUsersGuide.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard lo_density/*.eps) \
             $(wildcard flowchart/*.eps) \
             $(wildcard architecture/*.eps) \
             $(wildcard initial_models/*.eps) \
             $(wildcard rotation/*.eps) \
             $(wildcard radial_base/*.eps) \
             $(wildcard spherical_basestate/*.eps) \
             $(wildcard visualization/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = maestro_doc.bib \
                  introduction/introduction.tex \
                  flowchart/flowchart.tex \
                  preface/preface.tex \
                  architecture/architecture.tex \
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
                  analysis/analysis.tex \
                  thermo_notes/thermo_notes.tex \
                  visualization/visualization.tex


MaestroUsersGuide.dvi: MaestroUsersGuide.tex $(eps_source) $(section_sources) $(eps_source) 
	latex MaestroUsersGuide.tex < /dev/null
	bibtex MaestroUsersGuide
	latex MaestroUsersGuide.tex < /dev/null
	latex MaestroUsersGuide.tex < /dev/null
	latex MaestroUsersGuide.tex < /dev/null

MaestroUsersGuide.pdf: MaestroUsersGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex MaestroUsersGuide.tex < /dev/null
	bibtex MaestroUsersGuide
	pdflatex MaestroUsersGuide.tex < /dev/null
	pdflatex MaestroUsersGuide.tex < /dev/null
	pdflatex MaestroUsersGuide.tex < /dev/null

pdf:	MaestroUsersGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) MaestroUsersGuide.dvi 
	$(RM) MaestroUsersGuide.blg MaestroUsersGuide.log 
	$(RM) MaestroUsersGuide.toc MaestroUsersGuide.lof
	$(RM) MaestroUsersGuide.aux MaestroUsersGuide.ps MaestroUsersGuide.bbl
	$(RM) *~

.PHONY: clean
