# makefile for the element_sedimentation unit tests
#
# Builds ../tester from the solver modules in ../../code and the test sources in
# ../src against the MESA SDK libraries in $(MESA_DIR)/lib. Nothing is installed
# into MESA; MESA only has to be built (libconst, libmath, libutils).
# Object and module files are kept in this directory.

ifndef MESA_DIR
$(error MESA_DIR environment variable is not set)
endif

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

# only the MESA include directory is needed (no package-local public/make directories)
TEST_INCLUDES = -I$(MESA_INCLUDE_DIR)

#################################################################
#
# SOURCES

CODE_DIR = ../../code
SRC_DIR = ../src

MODULE_SRCS = advection_diffusion.f90 sedimentation_phase_diagram.f90 sedimentation_miscibility.f90
TEST_SRCS = test_support.f90 test_advection_diffusion.f90 test_phase_diagram.f90 test_miscibility.f90 tester.f90

SRCS = $(MODULE_SRCS) $(TEST_SRCS)
OBJS = $(patsubst %.f90,%.o,$(SRCS))

#################################################################
#
# LIBRARIES

LIBS_OTHER = utils math const
LOAD_OTHER = -L$(MESA_LIB_DIR) $(addprefix -l,$(LIBS_OTHER)) $(LOAD_CRMATH) $(LOAD_LAPACK) $(LOAD_BLAS)

#################################################################
#
# TARGETS

TEST = $(TEST_DIR)/tester

all : $(TEST)

$(TEST) : $(OBJS)
ifneq ($(QUIET),)
	@echo LOADER $@
	@$(LOADER) $(FCopenmp) -o $@ $(OBJS) $(LOAD_OTHER) $(LD_FLAGS)
else
	$(LOADER) $(FCopenmp) -o $@ $(OBJS) $(LOAD_OTHER) $(LD_FLAGS)
endif

clean :
	-@rm -f *.o *.mod *.smod $(TEST)

#################################################################
#
# COMPILATION RULES

%.o : %.f90
ifneq ($(QUIET),)
	@echo TEST_COMPILE $<
	@$(TEST_COMPILE) $(FCfree) $<
else
	$(TEST_COMPILE) $(FCfree) $<
endif

#################################################################
#
# DEPENDENCIES

vpath %.f90 $(CODE_DIR) $(SRC_DIR)

# module dependencies (compile order); listed explicitly instead of using makedepf90
sedimentation_miscibility.o : sedimentation_phase_diagram.o
test_advection_diffusion.o : advection_diffusion.o test_support.o
test_phase_diagram.o : sedimentation_phase_diagram.o test_support.o
test_miscibility.o : sedimentation_miscibility.o sedimentation_phase_diagram.o test_support.o
tester.o : test_support.o test_advection_diffusion.o test_phase_diagram.o test_miscibility.o
