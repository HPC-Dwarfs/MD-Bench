#CONFIGURE BUILD SYSTEM
# Optional directory prefix for TARGET/BUILD_DIR (binary, objects, and generated .s files)
PREFIX_DIR ?=
override PREFIX_DIR := $(patsubst %/,%,$(strip $(PREFIX_DIR)))
TAG = $(OPT_TAG)-$(TOOL_TAG)-$(DATA_TYPE)$(TAG_SUFFIX)
TARGET     = $(if $(PREFIX_DIR),$(PREFIX_DIR)/,)MDBench-$(TAG)
BUILD_DIR  = $(if $(PREFIX_DIR),$(PREFIX_DIR)/,)build/build-$(TAG)
SRC_ROOT   = src
SRC_DIR    = $(SRC_ROOT)/$(OPT_SCHEME)
COMMON_DIR = $(SRC_ROOT)/common
MAKE_DIR   = ./make
Q         ?= @

#DO NOT EDIT BELOW
include config.mk
include $(MAKE_DIR)/include_$(TOOLCHAIN).mk
include $(MAKE_DIR)/include_LIKWID.mk
INCLUDES  += -I$(CURDIR)/$(SRC_DIR) -I$(CURDIR)/$(COMMON_DIR)
VPATH     = $(SRC_DIR) $(COMMON_DIR) $(CUDA_DIR)
ASM       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.s,$(wildcard $(SRC_DIR)/*.c))
OBJ       = $(filter-out $(BUILD_DIR)/main%, $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c)))
OBJ      += $(patsubst $(COMMON_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(COMMON_DIR)/*.c))
ifneq ($(filter $(strip $(TOOLCHAIN)), NVCC HIPCC),)
OBJ      += $(patsubst $(COMMON_DIR)/%.cu, $(BUILD_DIR)/%.o,$(wildcard $(COMMON_DIR)/*.cu))
OBJ      += $(patsubst $(SRC_DIR)/%.cu, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cu))
endif
SOURCES   =  $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*.c $(COMMON_DIR)/*.c $(COMMON_DIR)/*.h)
CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(OPTIONS) $(INCLUDES)
c := ,
clist = $(subst $(eval) ,$c,$(strip $1))

define CLANGD_TEMPLATE
CompileFlags:
  Add: [$(call clist,$(INCLUDES)), $(call clist,$(CPPFLAGS))]
  Compiler: clang
endef

ifneq ($(VARIANT),)
	.DEFAULT_GOAL := ${TARGET}-$(VARIANT)
    DEFINES += -DVARIANT=$(VARIANT)
endif

${TARGET}: $(BUILD_DIR) .clangd $(OBJ) $(SRC_DIR)/main.c
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${LINKER} $(CPPFLAGS) ${LFLAGS} -o $(TARGET) $(OBJ) $(SRC_DIR)/main.c $(LIBS)

${TARGET}-%: $(BUILD_DIR) $(OBJ) $(SRC_DIR)/main-%.c
	@echo "===>  LINKING  $(TARGET)-$* "
	$(Q)${LINKER} $(CPPFLAGS) ${LFLAGS} -o $(TARGET)-$* $(OBJ) $(SRC_DIR)/main-$*.c $(LIBS)

$(BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TOOLCHAIN).mk
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) $(CFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d

ifeq ($(strip $(TOOLCHAIN)),NVCC)
$(BUILD_DIR)/%.o:  %.cu
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $(OPT_CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) $(OPT_CFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d
endif

ifeq ($(strip $(TOOLCHAIN)),HIPCC)
$(BUILD_DIR)/%.o:  $(BUILD_DIR)/%.hip
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d
endif

$(BUILD_DIR)/%.hip: %.cu
	$(info ===>  GENERATE HIP  $@)
	$(Q)hipify-perl $< > $@

$(BUILD_DIR)/%.s:  %.c
	$(info ===>  GENERATE ASM  $@)
	$(Q)$(CC) -S $(ASFLAGS) $(CPPFLAGS) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o:  %.s
	$(info ===>  ASSEMBLE  $@)
	$(Q)$(AS) $< -o $@

.PHONY: clean distclean cleanall tags format info asm sass test

clean:
	$(info ===>  CLEAN)
	@rm -rf $(BUILD_DIR)

cleanall:
	$(info ===>  CLEAN ALL)
	@rm -rf build
	@rm -rf MDBench-*
	@rm -f tags .clangd

distclean: clean
	$(info ===>  DIST CLEAN)
	@rm -f $(TARGET)
	@rm -f tags .clangd

info:
	$(info $(CFLAGS))
	$(Q)$(CC) $(VERSION)

asm:  $(BUILD_DIR) $(ASM)

# GPU machine-code dump for the compiled kernels: NVIDIA SASS (cuobjdump) or
# AMD GCN/RDNA ISA (roc-obj, from the rocm-utils package).
ifeq ($(strip $(TOOLCHAIN)),NVCC)
sass: $(TARGET)
	$(info ===>  DUMP SASS  $(TARGET).sass)
	$(Q)cuobjdump --dump-sass $(TARGET) > $(TARGET).sass
else ifeq ($(strip $(TOOLCHAIN)),HIPCC)
sass: $(TARGET)
	$(info ===>  DUMP GCN ISA  $(TARGET).gcn)
	$(Q)roc-obj -o $(BUILD_DIR)/roc-obj $(TARGET)
	$(Q)cat $(BUILD_DIR)/roc-obj/*.s > $(TARGET).gcn
else
sass:
	$(error sass target requires TOOLCHAIN=NVCC or TOOLCHAIN=HIPCC)
endif

tags:
	$(info ===>  GENERATE  TAGS)
	$(Q)ctags -R

format:
	@for src in $(SOURCES) ; do \
		echo "Formatting $$src" ; \
		clang-format -i $$src ; \
	done
	@echo "Done"

TEST_BIN := mdbench-tests

test: $(TEST_BIN) $(TARGET)
	@echo "===>  RUNNING  $(TEST_BIN)"
	@./$(TEST_BIN)
	@echo "===>  RUNNING  sim_argon_regression on $(TARGET)"
	@bash tests/sim_argon_regression.sh ./$(TARGET)
	@echo "===>  RUNNING  sim_copper_fcc_regression on $(TARGET)"
	@bash tests/sim_copper_fcc_regression.sh ./$(TARGET)
	@echo "===>  RUNNING  regression_energy_lj on $(TARGET)"
	@bash tests/regression_energy_lj.sh ./$(TARGET)
	@echo "===>  RUNNING  regression_lj_table on $(TARGET)"
	@bash tests/regression_lj_table.sh ./$(TARGET)
	@echo "===>  RUNNING  test_double_cutoff on $(TARGET)"
	@bash tests/test_double_cutoff.sh ./$(TARGET)
	@echo "===>  RUNNING  regression_scheme_equiv (geometric)"
	@bash tests/regression_scheme_equiv.sh
	@echo "===>  RUNNING  regression_scheme_equiv (single)"
	@LJ_COMB_RULE=single bash tests/regression_scheme_equiv.sh
	@echo "===>  RUNNING  test_simd_vs_scalar"
	@bash tests/test_simd_vs_scalar.sh
	@echo "===>  RUNNING  test_half_neigh"
	@bash tests/test_half_neigh.sh
	@echo "===>  RUNNING  test_data_layout"
	@bash tests/test_data_layout.sh
	@echo "===>  RUNNING  test_mpi"
	@bash tests/test_mpi.sh
	@echo "===>  RUNNING  test_openmp"
	@bash tests/test_openmp.sh
	@echo "===>  RUNNING  test_displacement_reneigh"
	@bash tests/test_displacement_reneigh.sh
	@echo "===>  RUNNING  test_main_stub"
	@bash tests/test_main_stub.sh

TEST_COMMON_SRCS := $(COMMON_DIR)/parameter.c $(COMMON_DIR)/box.c $(COMMON_DIR)/thermo.c $(COMMON_DIR)/allocate.c $(COMMON_DIR)/util.c $(COMMON_DIR)/ljtable.c
TEST_CP_SRCS     := $(SRC_ROOT)/clusterpair/atom.c $(SRC_ROOT)/clusterpair/neighbor.c $(SRC_ROOT)/clusterpair/pbc.c $(SRC_ROOT)/clusterpair/integrate.c

$(TEST_BIN): tests/main.c tests/test_runner.h tests/test_parameter.c tests/test_atom.c tests/test_force.c tests/test_neighbor.c tests/test_integrate.c tests/test_box.c tests/test_thermo.c $(TEST_COMMON_SRCS) $(TEST_CP_SRCS)
	@echo "===>  BUILDING $(TEST_BIN)"
	$(Q)$(CC) -I$(CURDIR)/src/clusterpair -I$(CURDIR)/src/common $(CPPFLAGS) $(CFLAGS) -DCLUSTER_PAIR -DCLUSTERPAIR_KERNEL_AUTO \
		tests/main.c tests/test_parameter.c tests/test_atom.c tests/test_force.c tests/test_neighbor.c \
		tests/test_integrate.c tests/test_box.c tests/test_thermo.c \
		$(TEST_COMMON_SRCS) $(TEST_CP_SRCS) \
		-lm -o $(TEST_BIN)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

.clangd:
	$(file > .clangd,$(CLANGD_TEMPLATE))

-include $(OBJ:.o=.d)
