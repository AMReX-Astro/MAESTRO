f90sources += network.f90
f90sources += burner.f90

# network.f90 is created at build time for this network
network.f90:
	@echo "---------------------------------------------------------------------------"
	@echo "WRITING network.f90:"
	$(FPARALLEL)/extern/networks/general_null/write_network.py \
            -t $(FPARALLEL)/extern/networks/general_null/network.template \
            -s $(FPARALLEL)/extern/networks/general_null/simple.net \
            -o network.f90
	@echo "---------------------------------------------------------------------------"
	@echo " "


# remove network.f90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) network.f90
