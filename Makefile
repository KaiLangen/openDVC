#############################################################################
# Constants
#############################################################################

ENC_EXE    := enDVC
DEC_EXE    := deDVC
COL_EXE		 := colour
TST_EXE		 := tests

OUT_DIR    := bin

ENC_TARGET := $(OUT_DIR)/$(ENC_EXE)
DEC_TARGET := $(OUT_DIR)/$(DEC_EXE)
COL_TARGET := $(OUT_DIR)/$(COL_EXE)
TST_TARGET := $(OUT_DIR)/$(TST_EXE)

find_objs   = $(subst src/,obj/, \
              $(subst .cpp,.o, \
							$(foreach i,$1/src,$(wildcard $i/*.cpp)))) $(COM_OBJS)

COM_OBJS   := $(call find_objs,common)
ENC_OBJS   := $(call find_objs,encoder)
DEC_OBJS   := $(call find_objs,decoder)
COL_OBJS   := $(call find_objs,colourizer)
TST_OBJS   := $(call find_objs,tests)

# Flags
CXXFLAGS   := -Wall -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-variable
#CXXFLAGS   += -O3 
CXXFLAGS   += -g3 -std=c++11

# Utilities
SHELL      := /bin/sh
CXX        := g++
SED        := sed
RM         := rm -f
TAR        := /bin/tar
MAKE       := make
CD         := cd
CP         := cp -f


#############################################################################
# Rules
#############################################################################

ifneq "$(SUB_MAKE)" "yes"

.DEFAULT_GOAL := default

.PHONY: default
default: encoder decoder tests #colourizer bin
	@echo ""

.PHONY: bin
bin:
	@echo ""
	@echo "Building directory bin"
	@$(TAR) -xzvf bin/ldpca.tar.gz --directory bin

.PHONY: common
common:
	@echo ""
	@echo "Building directory common"
	@$(MAKE) --no-print-directory -s -C common

.PHONY: encoder
encoder: common
	@echo ""
	@echo "Building directory encoder"
	@$(MAKE) --no-print-directory -s -C encoder
	@$(CXX) $(ENC_OBJS) -o $(ENC_TARGET)
	@echo "  LD   $(ENC_EXE)"

.PHONY: decoder
decoder: common
	@echo ""
	@echo "Building directory decoder"
	@$(MAKE) --no-print-directory -s -C decoder
	@$(CXX) $(DEC_OBJS) -o $(DEC_TARGET)
	@echo "  LD   $(DEC_EXE)"

.PHONY: colourizer
colourizer: common
	@echo ""
	@echo "Building directory colourizer"
	@$(MAKE) --no-print-directory -s -C colourizer
	@$(CXX) $(COL_OBJS) -o $(COL_TARGET)
	@echo "  LD   $(COL_EXE)"

.PHONY: tests
tests: common
	@echo ""
	@echo "Building directory tests"
	@$(MAKE) --no-print-directory -s -C tests
	@$(CXX) $(TST_OBJS) -o $(TST_TARGET)
	@echo "  LD   $(TST_EXE)"

.PHONY: clean
clean:
	@$(RM) $(ENC_TARGET) $(DEC_TARGET) $(COL_TARGET) $(TST_TARGET) 
	@echo "Cleaning directory common"
	@$(MAKE) --no-print-directory -s -C common clean
	@echo "Cleaning directory encoder"
	@$(MAKE) --no-print-directory -s -C encoder clean
	@echo "Cleaning directory decoder"
	@$(MAKE) --no-print-directory -s -C decoder clean
	@echo "Cleaning directory colourizer"
	@$(MAKE) --no-print-directory -s -C colourizer clean
	@echo "Cleaning directory tests"
	@$(MAKE) --no-print-directory -s -C tests clean

.PHONY: clean_all
clean_all: clean
	@$(RM) bin/pattern* bin/wz.bin bin/rec.yuv bin/wz.y
	@$(RM) bin/jm/data.txt bin/jm/*.dat bin/jm/test.264
	@$(RM) -r bin/ldpca/

endif

