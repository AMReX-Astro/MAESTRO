EPStoPDF = epstopdf

ALL: MaestroUsersGuide.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard lo_density/*.eps) \
             $(wildcard flowchart/*.eps) \
             $(wildcard architecture/*.eps) \
             $(wildcard unit_tests/*.eps) \
             $(wildcard initial_models/*.eps) \
             $(wildcard rotation/*.eps) \
             $(wildcard radial_base/*.eps) \
             $(wildcard spherical_basestate/*.eps) \
             $(wildcard visualization/*.eps) \
             $(wildcard planar_invsq_base/*.eps) \
             $(wildcard pert_notes/*.eps) \
             $(wildcard SDC/*.eps) \
             $(wildcard test_suite/*.eps) \
             $(wildcard getting_started/*.eps) \
             $(wildcard mg/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = maestro_doc.bib \
                  introduction/introduction.tex \
                  flowchart/flowchart.tex \
                  preface/preface.tex \
                  getting_started/getting_started.tex \
                  architecture/architecture.tex \
                  runtime_parameters/runtime_parameters.tex \
                  volume_discrepancy/volume_discrepancy.tex \
                  eos_notes/eos_notes.tex \
                  networks/networks.tex \
                  eta_notes/eta.tex \
                  pert_notes/pert.tex \
                  Godunov_notes/Godunov.tex \
                  initial_models/initial_models.tex \
                  lo_density/lo_density.tex \
                  rotation/rotation.tex \
                  radial_base/radial_base.tex \
                  spherical_basestate/basestate.tex \
                  FAQ/faq.tex \
                  compilers/compilers.tex \
                  analysis/analysis.tex \
                  thermo_notes/thermo_notes.tex \
                  visualization/visualization.tex \
                  SDC/sdc.tex \
                  test_suite/test_suite.tex \
                  makefiles/makefiles.tex \
                  managing_jobs/managingjobs.tex \
                  mg/mg.tex


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
