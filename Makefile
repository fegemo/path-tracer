# points to NVIDIA's CUDA/OpenCL implementation
CUDASDKROOT=/usr/local/cuda

CC=gcc
CCFLAGS=-O3 -msse2 -mfpmath=sse -ftree-vectorize -funroll-loops -Wall \
	-I$(CUDASDKROOT)/include -L$(CUDASDKROOT)/lib64 -lm -lglut -lGL -lOpenCL

default: all

all: Makefile cpu gpu preprocessed_kernels

cpu: main-cpu.c displayfunc.c Makefile vec.h camera.h geom.h displayfunc.h simplernd.h geomfunc.h
	$(CC) -DSMALLPT_CPU -o cpu main-cpu.c displayfunc.c $(CCFLAGS)
	./cpu

gpu: main-gpu.c displayfunc.c Makefile vec.h camera.h geom.h displayfunc.h simplernd.h geomfunc.h
	$(CC) -DSMALLPT_GPU -o gpu main-gpu.c displayfunc.c $(CCFLAGS)
	./gpu

clean:
	rm -rf main-cpu main-gpu preprocessed_rendering_kernel.cl preprocessed_rendering_kernel_dl.cl

preprocessed_kernels:
	cpp <rendering_kernel.cl >preprocessed_rendering_kernel.cl
	cpp <rendering_kernel_dl.cl >preprocessed_rendering_kernel_dl.cl
