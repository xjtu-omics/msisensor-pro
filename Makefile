CXX ?= g++
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++11

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INSTALL_BIN ?= install -c -m 0755
INSTALL_DIR ?= install -d -m 0755

TARGET := msisensor-pro
MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
HTSLIB ?= $(abspath $(MAKEFILE_DIR)/vendor/htslib-1.24)
HTSLIB_LIB := $(HTSLIB)/libhts.a

# OS & OpenMP detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
    ifeq ($(LIBOMP_PREFIX),)
        ifneq ($(wildcard /opt/homebrew/opt/libomp),)
            LIBOMP_PREFIX := /opt/homebrew/opt/libomp
        else ifneq ($(wildcard /usr/local/opt/libomp),)
            LIBOMP_PREFIX := /usr/local/opt/libomp
        endif
    endif
    ifneq ($(LIBOMP_PREFIX),)
        OPENMP_CXXFLAGS ?= -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        OPENMP_LDFLAGS ?= -L$(LIBOMP_PREFIX)/lib -lomp
    else
        OPENMP_CXXFLAGS ?=
        OPENMP_LDFLAGS ?=
    endif
else
    OPENMP_CXXFLAGS ?= -fopenmp
    OPENMP_LDFLAGS ?= -fopenmp
endif

# Version extraction
VERSION ?=
ifeq ($(VERSION),)
    BASE_VERSION := $(strip $(shell cat $(MAKEFILE_DIR)/VERSION 2>/dev/null || echo "dev"))
    GIT_TAG := $(shell git describe --tags --exact-match 2>/dev/null)
    GIT_DESC := $(shell git describe --always --dirty 2>/dev/null)
    ifneq ($(GIT_TAG),)
        VERSION := $(GIT_TAG)
    else ifneq ($(GIT_DESC),)
        VERSION := $(BASE_VERSION)-$(GIT_DESC)
    else
        VERSION := $(BASE_VERSION)
    endif
endif

# Include HTSlib static flags if available
-include $(HTSLIB)/htslib_static.mk

CPPFLAGS += -I$(HTSLIB) $(OPENMP_CXXFLAGS) -DVERSION=\"$(VERSION)\"
LDFLAGS += $(OPENMP_LDFLAGS) $(HTSLIB_static_LDFLAGS)
LIBS += -pthread $(HTSLIB_static_LIBS)
ifeq ($(HTSLIB_static_LIBS),)
    LIBS += -lz -lm -lbz2 -llzma -lcurl
endif

SOURCE := $(sort $(wildcard cpp/*.cpp))
OBJS := $(SOURCE:.cpp=.o)
DEPS := $(SOURCE:.cpp=.d)

.PHONY: all binary install uninstall test check clean distclean htslib

all: $(TARGET)

htslib: $(HTSLIB_LIB)

$(HTSLIB)/configure:
	@if [ ! -f "$@" ]; then \
		echo "Generating HTSlib configure script..."; \
		cd $(HTSLIB) && (autoheader 2>/dev/null || true) && autoconf; \
	fi

$(HTSLIB)/Makefile: $(HTSLIB)/configure
	@if [ ! -f "$@" ]; then \
		echo "Configuring HTSlib..."; \
		cd $(HTSLIB) && ./configure; \
	fi

$(HTSLIB_LIB): $(HTSLIB)/Makefile
	@echo "Building HTSlib static library..."
	$(MAKE) -C $(HTSLIB) lib-static

%.o: %.cpp | $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -MP -c $< -o $@

$(TARGET): $(OBJS) $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(OBJS) $(HTSLIB_LIB) $(LDFLAGS) $(LIBS) -o $@

install: $(TARGET)
	$(INSTALL_DIR) $(DESTDIR)$(BINDIR)
	$(INSTALL_BIN) $(TARGET) $(DESTDIR)$(BINDIR)/$(TARGET)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET)

test check: $(TARGET)
	@echo "Running msisensor-pro smoke tests..."
	./$(TARGET) scan -d demo/data_input/reference/reference.fa -o .test_scan.list
	./$(TARGET) msi -d .test_scan.list -n demo/data_input/data4test/test_normal_sorted.bam -t demo/data_input/data4test/test_tumor_sorted.bam -g demo/data_input/reference/reference.fa -o .test_msi.out
	@rm -f .test_scan.list .test_msi.out*
	@echo "All tests passed successfully."

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET) .test_scan.list .test_msi.out*

distclean: clean
	rm -f binary/$(TARGET)
	@if [ -f "$(HTSLIB)/Makefile" ]; then $(MAKE) -C $(HTSLIB) clean; fi

-include $(DEPS)
