f90sources += network.f90
f90sources += burner.f90

# network.f90 is created at build time for this network
network.f90:   $(GENERAL_NET_INPUTS) $(MAESTRO_TOP_DIR)/Microphysics/networks/general_null/network.template
	@echo " "
	@echo "${bold}WRITING network.f90${normal}"
	$(MAESTRO_TOP_DIR)/Microphysics/networks/general_null/write_network.py \
            -t $(MAESTRO_TOP_DIR)/Microphysics/networks/general_null/network.template \
            -s $(GENERAL_NET_INPUTS) \
            -o network.f90
	@echo " "


# remove network.f90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) network.f90
