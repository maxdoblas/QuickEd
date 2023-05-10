###############################################################################
# Compile rules
###############################################################################
BUILD_RULES=bitpal_build \
            daligner_wave_build \
            edlib_build \
            parasail_build \
            scrooge_build
            
CLEAN_RULES=edlib_clean \
            parasail_clean

all: prepare $(BUILD_RULES)

prepare:
	-mkdir build

clean: $(CLEAN_RULES)
	-rm -rf build

###############################################################################
# Build rules
###############################################################################
bitpal_build:
	-$(MAKE) --directory=bitpal all

daligner_wave_build:
	-$(MAKE) --directory=daligner_wave all

edlib_build:
	-cd edlib; ./build.sh

parasail_build:
	-cd parasail; ./build.sh
	
scrooge_build:
	-$(MAKE) --directory=scrooge all

###############################################################################
# Clean rules
###############################################################################
edlib_clean:
	-cd edlib; ./clean.sh;

parasail_clean:
	-cd parasail; ./clean.sh;



