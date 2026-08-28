# LINUX is used in the sources as "not Windows" (OpenMP vs MSVC PPL).
# Keep defining it on Darwin so the OpenMP path is compiled.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
BINDIR = ../binaries

ifeq ($(UNAME_S),Darwin)
LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /opt/homebrew/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /usr/local/opt/libomp
endif
ifeq ($(UNAME_M),arm64)
LIBNAME = RenderIrregular_arm64.so
else
LIBNAME = RenderIrregular_x86_64.so
endif
# Apple clang rejects -fopenmp; Homebrew libomp needs the preprocessor flag.
CXXFLAGS = -c -std=c++11 -O3 -fPIC -D LINUX -Xpreprocessor -fopenmp -I$(LIBOMP)/include
LDFLAGS = -shared -L$(LIBOMP)/lib -lomp -Wl,-rpath,$(LIBOMP)/lib
else
LIBNAME = RenderIrregular.so
CXXFLAGS = -c -std=c++11 -O3 -fPIC -D LINUX -fopenmp
LDFLAGS = -shared -fopenmp
endif

ONAME = $(BINDIR)/$(LIBNAME)

RenderIrregular	:	Render.o ExtInterface.o ExtMath.o
			g++ $^ $(LDFLAGS) -o $(ONAME)
			rm -f Render.o ExtInterface.o ExtMath.o
Render.o	:	Render.cpp ExtMath.h Render.h
			g++ $(CXXFLAGS) Render.cpp
ExtInterface.o	:	ExtInterface.cpp Render.h
			g++ $(CXXFLAGS) ExtInterface.cpp
ExtMath.o	:	ExtMath.cpp
			g++ $(CXXFLAGS) ExtMath.cpp

.PHONY: install clean

install:	RenderIrregular

clean:
			rm -f *.o
