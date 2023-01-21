#!/usr/bin/make -f
# Makefile for MOD Convolution Reverb #
# ----------------------------------- #
# Created by falkTX
#

include dpf/Makefile.base.mk

all: plugins gen

# --------------------------------------------------------------

plugins:
	$(MAKE) all -C src

ifneq ($(CROSS_COMPILING),true)
gen: plugins dpf/utils/lv2_ttl_generator
	@$(CURDIR)/dpf/utils/generate-ttl.sh

dpf/utils/lv2_ttl_generator:
	$(MAKE) -C dpf/utils/lv2-ttl-generator
else
gen:
endif

# --------------------------------------------------------------

clean:
	$(MAKE) clean -C dpf/utils/lv2-ttl-generator
	$(MAKE) clean -C src
	rm -f 3rd-party/FFTConvolver/*.d 3rd-party/FFTConvolver/*.o

# --------------------------------------------------------------

.PHONY: plugins
