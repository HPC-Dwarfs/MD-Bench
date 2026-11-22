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
include $(MAKE_DIR)/options.mk
include $(MAKE_DIR)/include_$(TOOLCHAIN).mk
include $(MAKE_DIR)/include_LIKWID.mk
INCLUDES  += -I$(CURDIR)/$(SRC_DIR) -I$(CURDIR)/$(COMMON_DIR)
VPATH     = $(SRC_DIR) $(COMMON_DIR) $(CUDA_DIR)
ASM       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.s,$(wildcard $(SRC_DIR)/*.c))
OBJ       = $(filter-out $(BUILD_DIR)/main%, $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c)))
OBJ      += $(patsubst $(COMMON_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(COMMON_DIR)/*.c))
ifneq ($(filter $(strip $(TOOLCHAIN)), NVCC HIPCC),)
CU_OBJS   = $(patsubst $(COMMON_DIR)/%.cu, $(BUILD_DIR)/%.o,$(wildcard $(COMMON_DIR)/*.cu))
CU_OBJS  += $(patsubst $(SRC_DIR)/%.cu, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cu))
OBJ      += $(CU_OBJS)
SASS       = $(CU_OBJS:.o=.sass)
FORCE_SASS = $(foreach o,$(CU_OBJS),$(if $(findstring force_lj_cuda,$o),$(o:.o=-force.sass)))
endif
SOURCES   =  $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*.c $(COMMON_DIR)/*.c $(COMMON_DIR)/*.h)
FORMAT_SOURCES = $(shell find $(SRC_ROOT) tests -name '*.c' -o -name '*.h')
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
else
    # The binary (the kernel-variant rule below would otherwise be the first target)
    .DEFAULT_GOAL := $(TARGET)
endif

# Kernel variants (src/common/kernel_variant.h): a template file is compiled once per
# variant instead of once on its own, as <file>.<variant>.o. A variant is the LJ rule
# (single, geom, full) and options, joined by '-', each one adds the flags below.
VARIANT_FLAGS_single     = -DLJ_COMB_RULE=0
VARIANT_FLAGS_geom       = -DLJ_COMB_RULE=1
VARIANT_FLAGS_full       = -DLJ_COMB_RULE=2
VARIANT_FLAGS_accurate   = -DACCURATE_RECIPROCAL
VARIANT_FLAGS_noprefetch = -DKERNEL_NO_PREFETCH
VARIANT_FLAGS_excl       = -DKERNEL_EXCLUSIONS

RULES          = single geom full
RULES_EXCL     = $(foreach r,$(RULES),$(r) $(r)-excl)
RULES_ACCURATE = $(foreach r,$(RULES),$(r) $(r)-accurate)
# Rule x --accurate-recip x --no-prefetch x exclusions
RULES_SIMD     = $(foreach e,$(RULES_ACCURATE),$(foreach p,$(e) $(e)-noprefetch,$(p) $(p)-excl))

# The template files of the scheme and their variants (none: the file is not used)
ifeq ($(strip $(OPT_SCHEME)),clusterpair)
TEMPLATE_FILES          = neighbor_build force_lj force_lj_table force_lj_4xn force_lj_2xnn force_lj_2xn
VARIANTS_neighbor_build = $(RULES)
VARIANTS_force_lj       = $(RULES)
VARIANTS_force_lj_table = $(RULES)
ifneq ($(strip $(USE_REFERENCE_KERNEL)),true)
VARIANTS_force_lj_4xn   = $(RULES_SIMD)
VARIANTS_force_lj_2xnn  = $(RULES_SIMD)
VARIANTS_force_lj_2xn   = $(RULES_SIMD)
endif
else
TEMPLATE_FILES          = neighbor_build force_lj force_lj_table force_eam force_lj_simd
VARIANTS_neighbor_build = $(RULES)
VARIANTS_force_lj       = $(RULES_EXCL)
VARIANTS_force_lj_table = $(RULES)
VARIANTS_force_eam      = $(RULES)
ifeq ($(strip $(USE_SIMD_KERNEL)),true)
VARIANTS_force_lj_simd  = $(RULES_ACCURATE)
endif
endif

VARIANT_OBJ = $(foreach f,$(TEMPLATE_FILES),$(foreach v,$(VARIANTS_$(f)),$(BUILD_DIR)/$(f).$(v).o))
OBJ        := $(filter-out $(TEMPLATE_FILES:%=$(BUILD_DIR)/%.o),$(OBJ)) $(VARIANT_OBJ)
# Flags of a variant: $(call variant_flags,<file>.<variant>)
variant_flags = $(foreach o,$(subst -, ,$(subst .,,$(suffix $(1)))),$(VARIANT_FLAGS_$(o)))

# Extra flags of the toolchain for the SIMD force kernels (KERNEL_CFLAGS)
$(filter $(addprefix $(BUILD_DIR)/,$(addsuffix .%,force_lj_4xn force_lj_2xnn force_lj_2xn)),$(VARIANT_OBJ)): CFLAGS += $(KERNEL_CFLAGS)

.SECONDEXPANSION:
$(VARIANT_OBJ): $(BUILD_DIR)/%.o: $$(basename $$*).c $(MAKE_DIR)/include_$(TOOLCHAIN).mk $(CONFIG_STAMP)
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $(call variant_flags,$*) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) $(CFLAGS) $(call variant_flags,$*) -MT $@ -MM $< > $(@:.o=.d)

# Rebuild everything when the compiler or any flag changes (options outside the binary tag)
CONFIG_STAMP  = $(BUILD_DIR)/config.stamp
CONFIG_STRING = $(CC) $(LINKER) $(CPPFLAGS) $(CFLAGS) $(KERNEL_CFLAGS) $(OPT_CFLAGS) $(ASFLAGS) $(LFLAGS) $(LIBS)
ifeq ($(filter clean distclean cleanall format tags info print-target,$(MAKECMDGOALS)),)
ifneq ($(strip $(CONFIG_STRING)),$(strip $(file <$(CONFIG_STAMP))))
$(shell mkdir -p $(BUILD_DIR))
$(file >$(CONFIG_STAMP),$(CONFIG_STRING))
endif
endif

${TARGET}: $(BUILD_DIR) .clangd $(OBJ) $(SRC_DIR)/main.c $(CONFIG_STAMP)
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${LINKER} $(CPPFLAGS) ${LFLAGS} -o $(TARGET) $(OBJ) $(SRC_DIR)/main.c $(LIBS)

${TARGET}-%: $(BUILD_DIR) $(OBJ) $(SRC_DIR)/main-%.c $(CONFIG_STAMP)
	@echo "===>  LINKING  $(TARGET)-$* "
	$(Q)${LINKER} $(CPPFLAGS) ${LFLAGS} -o $(TARGET)-$* $(OBJ) $(SRC_DIR)/main-$*.c $(LIBS)

$(BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TOOLCHAIN).mk $(CONFIG_STAMP)
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) $(CFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d

ifeq ($(strip $(TOOLCHAIN)),NVCC)
$(BUILD_DIR)/%.o:  %.cu $(CONFIG_STAMP)
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $(OPT_CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) $(OPT_CFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d
endif

ifeq ($(strip $(TOOLCHAIN)),HIPCC)
$(BUILD_DIR)/%.o:  $(BUILD_DIR)/%.hip $(CONFIG_STAMP)
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $(HIP_DEVICE_CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -MT $@ -MM  $< > $(BUILD_DIR)/$*.d
endif

$(BUILD_DIR)/%.hip: %.cu
	$(info ===>  GENERATE HIP  $@)
	$(Q)hipify-perl $< > $@

$(BUILD_DIR)/%.s:  %.c $(CONFIG_STAMP)
	$(info ===>  GENERATE ASM  $@)
	$(Q)$(CC) -S $(ASFLAGS) $(CPPFLAGS) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o:  %.s
	$(info ===>  ASSEMBLE  $@)
	$(Q)$(AS) $< -o $@

.PHONY: clean distclean cleanall tags format info print-target asm sass test mdbench-tests

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
	@rm -f $(TARGET) $(TEST_BIN)
	@rm -f tags .clangd

info:
	$(info $(CFLAGS))
	$(Q)$(CC) $(VERSION)

# Path of the binary of a configuration, for scripts: make -s <options> print-target
print-target:
	@echo $(abspath $(TARGET))

asm:  $(BUILD_DIR) $(ASM)

# GPU machine-code dump for the compiled kernels, one file per .cu source
# (mirrors `asm`): NVIDIA SASS (cuobjdump) or AMD GCN/RDNA ISA (roc-obj, from
# the rocm-utils package). `sass-force` additionally filters each dump down
# to just the LJ force kernels (computeForceLJCuda*), for a quicker look at
# the hot kernels without wading through integrate/pbc/prune/pack kernels.
ifeq ($(strip $(TOOLCHAIN)),NVCC)
$(BUILD_DIR)/%.sass: $(BUILD_DIR)/%.o
	$(info ===>  DUMP SASS  $@)
	$(Q)cuobjdump --dump-sass $< > $@

$(BUILD_DIR)/%-force.sass: $(BUILD_DIR)/%.o
	$(info ===>  DUMP SASS (force kernels)  $@)
	$(Q)fn=$$(cuobjdump --dump-sass $< 2>/dev/null | grep 'Function :' | grep -i ForceLJ | awk '{print $$3}' | sort -u | paste -sd,); \
	if [ -n "$$fn" ]; then cuobjdump --function $$fn --dump-sass $< > $@; \
	else echo "warning: no ForceLJ kernel found in $<" >&2; : > $@; fi
else ifeq ($(strip $(TOOLCHAIN)),HIPCC)
$(BUILD_DIR)/%.sass: $(BUILD_DIR)/%.o
	$(info ===>  DUMP GCN ISA  $@)
	$(Q)rm -rf $@.tmp
	$(Q)roc-obj -d -o $@.tmp $< >/dev/null 2>&1; \
	cat $@.tmp/*.s > $@ 2>/dev/null || : > $@
	$(Q)rm -rf $@.tmp

$(BUILD_DIR)/%-force.sass: $(BUILD_DIR)/%.sass
	$(info ===>  FILTER GCN ISA (force kernels)  $@)
	$(Q)awk '/^[0-9a-f]+ </{keep=($$0 ~ /ForceLJ/)} keep' $< > $@
else
sass sass-force:
	$(error sass target requires TOOLCHAIN=NVCC or TOOLCHAIN=HIPCC)
endif

sass: $(BUILD_DIR) $(SASS)

sass-force: $(BUILD_DIR) $(FORCE_SASS)

tags:
	$(info ===>  GENERATE  TAGS)
	$(Q)ctags -R

format:
	@for src in $(FORMAT_SOURCES) ; do \
		echo "Formatting $$src" ; \
		clang-format -i $$src ; \
	done
	@echo "Done"

# Unit tests, tagged like the benchmark binary so every configuration gets its own
TEST_BIN = $(TARGET)-tests

mdbench-tests: $(TEST_BIN)

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
	@LJ_COMB=single bash tests/regression_scheme_equiv.sh
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
	@echo "===>  RUNNING  test_coulomb_exclusions"
	@bash tests/test_coulomb_exclusions.sh
	@echo "===>  RUNNING  test_mpi_gmxbench"
	@bash tests/test_mpi_gmxbench.sh
	@echo "===>  RUNNING  test_halfshell"
	@bash tests/test_halfshell.sh

TEST_COMMON_SRCS := $(COMMON_DIR)/parameter.c $(COMMON_DIR)/box.c $(COMMON_DIR)/thermo.c $(COMMON_DIR)/allocate.c $(COMMON_DIR)/util.c $(COMMON_DIR)/ljtable.c $(COMMON_DIR)/device.c
TEST_CP_SRCS     := $(SRC_ROOT)/clusterpair/atom.c $(SRC_ROOT)/clusterpair/neighbor.c $(SRC_ROOT)/clusterpair/pbc.c $(SRC_ROOT)/clusterpair/integrate.c

TEST_HEADERS     := $(wildcard $(SRC_ROOT)/clusterpair/*.h $(COMMON_DIR)/*.h $(COMMON_DIR)/simd/*.h)

# The cluster-pair neighbor-list build variants (one per LJ rule) for the unit tests
TEST_BUILD_DIR    = $(BUILD_DIR)/tests
TEST_NEIGHBOR_OBJ = $(RULES:%=$(TEST_BUILD_DIR)/neighbor_build.%.o)
$(TEST_NEIGHBOR_OBJ): $(TEST_BUILD_DIR)/neighbor_build.%.o: $(SRC_ROOT)/clusterpair/neighbor_build.c $(TEST_HEADERS) $(CONFIG_STAMP)
	@mkdir -p $(TEST_BUILD_DIR)
	$(Q)$(CC) -c -I$(CURDIR)/src/clusterpair -I$(CURDIR)/src/common $(CPPFLAGS) $(CFLAGS) -DCLUSTER_PAIR -DCLUSTERPAIR_KERNEL_AUTO \
		$(VARIANT_FLAGS_$*) $< -o $@

$(TEST_BIN): $(CONFIG_STAMP) $(TEST_NEIGHBOR_OBJ) tests/main.c tests/test_runner.h tests/test_parameter.c tests/test_atom.c tests/test_force.c tests/test_neighbor.c tests/test_integrate.c tests/test_box.c tests/test_thermo.c $(TEST_COMMON_SRCS) $(TEST_CP_SRCS) $(TEST_HEADERS)
	@mkdir -p $(dir $(TEST_BIN))
	@echo "===>  BUILDING $(TEST_BIN)"
	$(Q)$(CC) -I$(CURDIR)/src/clusterpair -I$(CURDIR)/src/common $(CPPFLAGS) $(CFLAGS) -DCLUSTER_PAIR -DCLUSTERPAIR_KERNEL_AUTO \
		tests/main.c tests/test_parameter.c tests/test_atom.c tests/test_force.c tests/test_neighbor.c \
		tests/test_integrate.c tests/test_box.c tests/test_thermo.c \
		$(TEST_COMMON_SRCS) $(TEST_CP_SRCS) $(TEST_NEIGHBOR_OBJ) \
		-lm -o $(TEST_BIN)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

.clangd:
	$(file > .clangd,$(CLANGD_TEMPLATE))

-include $(OBJ:.o=.d)
