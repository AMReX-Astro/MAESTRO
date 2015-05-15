# This is the makefile for testing the mesa net library

#Already defined in bashrc
#MESA_DIR = ../../.. 

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

#LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnet
LOAD_OTHER = -L$(MESA_LIB_DIR) -lnet -leos -lrates -lchem $(LOAD_MESA_NUMERICS)

#################################################################

# STEP 3: build

TEST_STD_OBJS = test_net_support.o mod_test_net.o
TEST = tester
TEST_OBJS = $(TEST_STD_OBJS) test_net.o

all : $(TEST)

DEBUG_FLAGS = -g

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)
	
#################################################################
   
%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfixed) $<

clean:
	-@rm -f *.o *.mod $(TEST_DIR)/$(TEST) 
