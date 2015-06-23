# This is the makefile for testing the mesa net library

#MESA_DIR = ../../..

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

ifeq ($(NET_OFFLOAD),YES)
CC += -Doffload
FC += -Doffload
endif

#################################################################

# STEP 2: info about the libraries

APROX21_FXT_FLAG =
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnet

#APROX21_FXT_FLAG = -DAPROX21
#LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnet -L$(MESA_DIR)/net/aprox21_fxt -la21fxt

#######

LOAD_OTHER = -L$(MESA_LIB_DIR) -leos -lrates -lchem $(LOAD_MESA_NUMERICS)

#################################################################

# STEP 3: build

TEST_BASE_OBJS = test_net_do_one.o test_net_support.o 
TEST_STD_OBJS = $(TEST_BASE_OBJS) mod_test_net.o

TEST = tester
TEST_OBJS = $(TEST_STD_OBJS) test_net.o

#PLOT = plotter
#PLOT_OBJS = test_net_support.o test_burn.o plot_net.o

all : $(TEST)

DEBUG_FLAGS = -g

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

	
#################################################################
   
TEST_COMPILE += $(APROX21_FXT_FLAG)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod $(TEST_DIR)/$(TEST) 
