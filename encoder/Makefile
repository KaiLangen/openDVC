
SUB_MAKE  := yes

-include ../Makefile

#############################################################################
# Constants
#############################################################################

# Source directories
SRC_DIR   := src
# Include directories
INC_DIR   := inc ../common/inc

# Source files
SRC_FILES := $(foreach i,$(SRC_DIR),$(wildcard $i/*.cpp))
# Object files
OBJ_FILES := $(subst .cpp,.o,$(filter %.cpp,$(SRC_FILES)))
OBJ_FILES := $(addprefix obj/,$(notdir $(OBJ_FILES)))
# Dependency files
DEP_FILES := $(subst .o,.d,$(OBJ_FILES))
DEP_FILES := $(addprefix dep/,$(notdir $(DEP_FILES)))

# Flags
CXXFLAGS  += $(addprefix -I,$(INC_DIR)) -DENCODER

VPATH		  := $(SRC_DIR)

ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEP_FILES)
endif

#############################################################################
# Rules
#############################################################################

.DEFAULT_GOAL := default

.PHONY: default
default: $(OBJ_FILES)

.PHONY: clean
clean:
	@$(RM) $(TARGET)
	@$(RM) $(OBJ_FILES)
	@$(RM) $(DEP_FILES)

obj/%.o: %.cpp
	@echo "  CXX  $(notdir $<)"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

dep/%.d: %.cpp
	@$(CXX) $(CXXFLAGS) -MM $< > $@.$$$$; \
	$(SED) 's,\($(notdir $*)\)\.o[ :]*,obj/\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$
