###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=../../../..
FOLDER_TOOL=../..
FOLDER_BUILD=$(FOLDER_TOOL)/build

###############################################################################
# Modules
###############################################################################
MODULES=benchmark_bitpal \
        benchmark_daligner \
        benchmark_diffutils \
        benchmark_edlib \
        benchmark_lv89 \
        benchmark_parasail \
        benchmark_scrooge_bindings \
        benchmark_scrooge \
        benchmark_seqan_bindings \
        benchmark_seqan
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

LI_FLAGS=-I$(FOLDER_ROOT) -I$(FOLDER_TOOL)

###############################################################################
# Rules
###############################################################################
all: $(OBJS) $(FOLDER_BUILD)/benchmark_seqan_bindings.o

$(FOLDER_BUILD)/benchmark_diffutils.o : benchmark_diffutils.c
	$(CC) -Wno-maybe-uninitialized $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@

$(FOLDER_BUILD)/benchmark_scrooge_bindings.o : benchmark_scrooge_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) $(LI_FLAGS) -I$(FOLDER_TOOL)/external/ -c $< -o $@

$(FOLDER_BUILD)/benchmark_seqan_bindings.o : benchmark_seqan_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) $(LI_FLAGS) -I$(FOLDER_TOOL)/external/seqan/include/ -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@
	
