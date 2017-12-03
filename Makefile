# points to NVIDIA's CUDA/OpenCL implementation
CUDASDKROOT=/usr/local/cuda

CC=gcc
CCFLAGS=-O3 -msse2 -mfpmath=sse -ftree-vectorize -funroll-loops -Wall \
	-I$(CUDASDKROOT)/include -L$(CUDASDKROOT)/lib64 -lm -lglut -lGL -lOpenCL

default: all

all: Makefile cpu gpu preprocessed_kernels

cpu: main-cpu.c displayfunc.c camera.c Makefile vec.h camera.h geom.h displayfunc.h simplernd.h geomfunc.h
	$(CC) -DSMALLPT_CPU -o path-tracer-cpu main-cpu.c displayfunc.c camera.c $(CCFLAGS)

gpu: main-gpu.c displayfunc.c camera.c Makefile vec.h camera.h geom.h displayfunc.h simplernd.h geomfunc.h
	$(CC) -DSMALLPT_GPU -o path-tracer-gpu main-gpu.c displayfunc.c camera.c $(CCFLAGS)

run-cpu: cpu
	./path-tracer-cpu

run-gpu: gpu
	./path-tracer-gpu

clean:
	rm -rf path-tracer-cpu path-tracer-gpu preprocessed_path_tracing.cl preprocessed_ray_tracing.cl

preprocessed_kernels:
	cpp <path-tracing.cl >preprocessed_path_tracing.cl
	cpp <ray-tracing.cl >preprocessed_ray_tracing.cl
