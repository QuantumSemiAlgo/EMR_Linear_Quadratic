# =================== Paths ===================
# Auto-detect or use environment variables
# Set PETSC_DIR and SLEPC_DIR in your environment, or they will be auto-detected
PETSC_DIR   ?= $(shell if [ -d "${HOME}/petsc" ]; then echo "${HOME}/petsc"; \
                       elif [ -d "/usr/lib/petsc" ]; then echo "/usr/lib/petsc"; \
                       elif [ -d "/opt/petsc" ]; then echo "/opt/petsc"; \
                       else echo "PETSC_NOT_FOUND"; fi)
SLEPC_DIR   ?= $(shell if [ -d "${HOME}/slepc" ]; then echo "${HOME}/slepc"; \
                       elif [ -d "/usr/lib/slepc" ]; then echo "/usr/lib/slepc"; \
                       elif [ -d "/opt/slepc" ]; then echo "/opt/slepc"; \
                       else echo "SLEPC_NOT_FOUND"; fi)

# =================== Compiler ===================
MPICXX       ?= mpicxx
CXX          := $(MPICXX)
ASAN_FLAGS   ?= -g -O1 -fsanitize=address -fno-omit-frame-pointer
# include paths
INC_FLAGS    = -I$(PETSC_DIR)/include \
               -I$(PETSC_DIR)/include/petsc \
               -I$(SLEPC_DIR)/include
# combine
CXXFLAGS     = $(INC_FLAGS) $(ASAN_FLAGS) -std=c++14

# =================== Files ===================
FILES        = readin_1 input_reader_1 main_1 utilrect_1 make_global_1 apply_bc_1 petsc_diagonalizer_1 solution_clean_1 convertnode_1 mesh_generator_1
OBJFILES     = $(addsuffix .o, $(FILES))
BINARY_NAME  = Laplace

# =================== Linker Flags ===================
LDFLAGS      = -L$(PETSC_DIR)/lib \
               -L$(SLEPC_DIR)/lib \
               -lslepc -lpetsc \
               $(ASAN_FLAGS)

# =================== Targets ===================
all: check-deps $(BINARY_NAME)

check-deps:
	@if [ "$(PETSC_DIR)" = "PETSC_NOT_FOUND" ]; then \
		echo "ERROR: PETSc not found. Set PETSC_DIR environment variable."; \
		exit 1; \
	fi
	@if [ "$(SLEPC_DIR)" = "SLEPC_NOT_FOUND" ]; then \
		echo "ERROR: SLEPc not found. Set SLEPC_DIR environment variable."; \
		exit 1; \
	fi
	@echo "Using PETSC_DIR=$(PETSC_DIR)"
	@echo "Using SLEPC_DIR=$(SLEPC_DIR)"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BINARY_NAME): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $(OBJFILES) $(LDFLAGS) -o $@

.PHONY: clean remake print-vars check-deps

clean:
	$(RM) $(OBJFILES) $(BINARY_NAME) *~

remake: clean all

print-vars:
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "LDFLAGS   = $(LDFLAGS)"
	@echo "CXX       = $(CXX)"
