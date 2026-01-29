# =================== Paths ===================
PETSC_DIR   ?= /Users/aryaer/Desktop/petsc-install
SLEPC_DIR   ?= /Users/aryaer/Desktop/slepc-install

# =================== Compiler ===================
# =================== Compiler ===================
CXX          ?= mpicxx
ASAN_FLAGS   ?= -g -O1 -fsanitize=address -fno-omit-frame-pointer
# include paths
INC_FLAGS    = -I$(PETSC_DIR)/include \
               -I$(PETSC_DIR)/include/petsc \
               -I$(SLEPC_DIR)/include
# combine
CXXFLAGS     = $(INC_FLAGS) $(ASAN_FLAGS)

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
all: $(BINARY_NAME)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BINARY_NAME): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $(OBJFILES) $(LDFLAGS) -o $@

.PHONY: clean remake print-vars

clean:
	$(RM) $(OBJFILES) $(BINARY_NAME) *~

remake: clean all

print-vars:
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "LDFLAGS   = $(LDFLAGS)"
	@echo "CXX       = $(CXX)"
