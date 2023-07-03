#!/usr/bin/make -f
# Makefile for DISTRHO Plugins #
# ---------------------------- #
# Created by falkTX
#

# --------------------------------------------------------------
# Project name, used for binaries

NAME = MOD-ConvolutionLoaderMono

# --------------------------------------------------------------
# Files to build

FILES_DSP = \
	OneKnobPlugin.cpp \
	3rd-party.cpp \
	../3rd-party/FFTConvolver/AudioFFT.cpp \
	../3rd-party/FFTConvolver/FFTConvolver.cpp \
	../3rd-party/FFTConvolver/TwoStageFFTConvolver.cpp \
	../3rd-party/FFTConvolver/Utilities.cpp

# --------------------------------------------------------------
# Do some magic

include ../dpf/Makefile.plugins.mk

BUILD_CXX_FLAGS += -DAUDIOFFT_PFFFT
BUILD_CXX_FLAGS += -DR8B_FASTTIMING
BUILD_CXX_FLAGS += -DR8B_R8B_PFFFT
BUILD_CXX_FLAGS += -I../3rd-party
BUILD_CXX_FLAGS += -I../3rd-party/FFTConvolver
BUILD_CXX_FLAGS += -I../3rd-party/r8brain
BUILD_CXX_FLAGS += -pthread
LINK_FLAGS      += $(SHARED_MEMORY_LIBS)

# BUILD_CXX_FLAGS += -DAUDIOFFT_FFTW3 $(shell $(PKG_CONFIG) --cflags fftw3f)
# LINK_FLAGS      += $(shell $(PKG_CONFIG) --libs fftw3f)

# --------------------------------------------------------------
# Enable all relevant plugin types

all: lv2_sep

# --------------------------------------------------------------
