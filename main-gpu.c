#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <CL/cl.h>

#include "camera.h"
#include "vec.h"
#include "geom.h"
#include "simplernd.h"
#include "displayfunc.h"
#include "time-utils.h"

/* Options */
static int useGPU = 1;
static int forceWorkSize = 0;

/* OpenCL variables */
static cl_context context;
static cl_mem colorBuffer;
static cl_mem pixelBuffer;
static cl_mem seedBuffer;
static cl_mem objectBuffer;
static cl_mem cameraBuffer;
static cl_mem debugBuffer;
static cl_command_queue commandQueue;
static cl_program program;
static cl_kernel kernel;
static unsigned int workGroupSize = 1;
static char *kernelFileName = "path-tracing.cl";
double startRenderingTime;

static vec *colors;
static unsigned int *seeds;
Camera camera;
int currentSample = 0;
Object *objects;
unsigned int objectCount;
unsigned int lightCount;
int *debug;


extern char captionLine1[];
extern char captionLine2[];


void updateRenderingStatistics(double, int);

///
/// Frees the color, pixel and random seeds buffer from vram and ram.
///
static void freeBuffers() {
	cl_int status = clReleaseMemObject(colorBuffer);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to release OpenCL color buffer: %d\n", status);
		exit(-1);
    }

	status = clReleaseMemObject(pixelBuffer);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to release OpenCL pixel buffer: %d\n", status);
		exit(-1);
    }

	status = clReleaseMemObject(seedBuffer);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to release OpenCL seed buffer: %d\n", status);
		exit(-1);
    }

	free(seeds);
	free(colors);
	free(pixels);
}

///
/// Allocates the output buffers necessary to show the partial results using glut.
///
static void allocateOutputBuffers() {
	const int pixelCount = width * height;
	int i;
	colors = (vec *)malloc(sizeof(vec) * pixelCount);

	// generates 2 random numbers per pixel (not sure what it is used for)
	seeds = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount * 2);
	for (i = 0; i < pixelCount * 2; i++) {
		seeds[i] = rand();
		if (seeds[i] < 2)
			seeds[i] = 2;
	}

	// array of indices of pixels (not sure what for)
	pixels = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount);
	for (i = 0; i < pixelCount; ++i) {
		pixels[i] = i;
	}

    debug = (int*)malloc(10*sizeof(int));

	// creates the color buffer, from which the opencl program can write (results) and read (what for?)
	cl_int status;
	cl_uint sizeBytes = sizeof(vec) * width * height;
    colorBuffer = clCreateBuffer(
            context,
            CL_MEM_READ_WRITE,
            sizeBytes,
            NULL,
            &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL output buffer: %d\n", status);
		exit(-1);
    }

    // creates a buffer of pixel indices, which the opencl can only read (what for?)
	sizeBytes = sizeof(unsigned int) * width * height;
    pixelBuffer = clCreateBuffer(
            context,
            CL_MEM_WRITE_ONLY,
            sizeBytes,
            NULL,
            &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL pixel buffer: %d\n", status);
		exit(-1);
    }

    sizeBytes = sizeof(int) * 10;
    debugBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeBytes, NULL, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL debugs buffer: %d\n", status);
		exit(-1);
    }

    // creates a buffer with random seeds (2 x pixel), read/write (what for?)
	sizeBytes = sizeof(unsigned int) * width * height * 2;
	seedBuffer = clCreateBuffer(
            context,
            CL_MEM_READ_WRITE,
            sizeBytes,
            NULL,
            &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL seed buffer: %d\n", status);
		exit(-1);
    }

    // writes the random seeds buffer to the vram
	status = clEnqueueWriteBuffer(
			commandQueue,
			seedBuffer,
			CL_TRUE,
			0,
			sizeBytes,
			seeds,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to write the OpenCL seeds buffer: %d\n", status);
		exit(-1);
	}
}

///
/// Reads a file, by its file name, line by line and returns an array of null-terminated lines.
///
static char *readFileByLine(const char *fileName) {
	FILE *file = fopen(fileName, "r");
	if (!file) {
		fprintf(stderr, "Failed to open file '%s'\n", fileName);
		exit(-1);
	}

	if (fseek(file, 0, SEEK_END)) {
		fprintf(stderr, "Failed to seek file '%s'\n", fileName);
		exit(-1);
	}

	long size = ftell(file);
	if (size == 0) {
		fprintf(stderr, "Failed to check position on file '%s'\n", fileName);
		exit(-1);
	}

	rewind(file);

	char *src = (char *)malloc(sizeof(char) * size + 1);
	if (!src) {
		fprintf(stderr, "Failed to allocate memory for file '%s'\n", fileName);
		exit(-1);
	}

	fprintf(stderr, "Reading file '%s' (size %ld bytes)\n", fileName, size);
	size_t res = fread(src, 1, sizeof(char) * size, file);
	if (res != sizeof(char) * size) {
		fprintf(stderr, "Failed to read file '%s' (read %ld)\n", fileName, res);
		exit(-1);
	}
	src[size] = '\0'; /* NULL terminated */

	fclose(file);

	return src;

}

///
/// Sets up platorm, devices, context, command queue, buffers and kernel.
///
static void setUpOpenCL() {
	// selecting the platform
    // ----------------------

    cl_uint numPlatforms;
	cl_platform_id platform = NULL;

	// num_entries (amount of platforms we want), platforms (their ids), num_platorms (amount of available)
	// in this case, we want to know how many we have only (hence [0, NULL, &numPlatforms])
	// so we can later fetch the details of the first one
	cl_int status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to get OpenCL platforms\n");
		exit(-1);
	}

	if (numPlatforms > 0) {
		cl_platform_id *platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id) * numPlatforms);
		// now we query again to fetch the available platforms' ids...
		status = clGetPlatformIDs(numPlatforms, platforms, NULL);
			if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL platform IDs\n");
			exit(-1);
		}

		unsigned int i;
		for (i = 0; i < numPlatforms; ++i) {
			char pbuf[100];
			// now we query details for each platform, in particular, its vendor (eg, NVIDIA)
			status = clGetPlatformInfo(platforms[i],
					CL_PLATFORM_VENDOR,
					sizeof(pbuf),
					pbuf,
					NULL);

			fprintf(stderr, "OpenCL Platform %d: %s\n", i, pbuf);
		}

		// we've just listed (printed) all of the available platforms...
		// now we simply get the first platform that was returned and that's it
		platform = platforms[0];
		free(platforms);
	}


	// selecting the device
	// --------------------

	cl_device_id devices[32];
	cl_uint deviceCount;
	// we query which available devices this platform has...
	// platform (id), device_type (cpu, gpu etc.), num_entries (max devices), devices (array of ids),
	//      num_devices (actual number of devices found on this platform)
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 32, devices, &deviceCount);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to get OpenCL device IDs\n");
		exit(-1);
	}

	int deviceFound = 0;
	cl_device_id selectedDevice;
	unsigned int i;
	for (i = 0; i < deviceCount; ++i) {
		cl_device_type type = 0;

		// query which type this device has and puts it on type...
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_TYPE,
				sizeof(cl_device_type),
				&type,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		char *stype;
		// generates a string with the device type so we can printf it
		switch (type) {
			case CL_DEVICE_TYPE_ALL:
				stype = "TYPE_ALL";
				break;
			case CL_DEVICE_TYPE_DEFAULT:
				stype = "TYPE_DEFAULT";
				break;
			case CL_DEVICE_TYPE_CPU:
				stype = "TYPE_CPU";
			    // we found a CPU, so we can use it and proceed with this program (if we're !useGPU)
				if (!useGPU && !deviceFound) {
					selectedDevice = devices[i];
					deviceFound = 1;
				}
				break;
			case CL_DEVICE_TYPE_GPU:
			    // we found a GPU, so we can use it and proceed with this program
				stype = "TYPE_GPU";
				if (useGPU && !deviceFound) {
					selectedDevice = devices[i];
					deviceFound = 1;
				}
				break;
			default:
				stype = "TYPE_UNKNOWN";
				break;
		}
		fprintf(stderr, "OpenCL Device %d: Type = %s\n", i, stype);

		// now we query the name of this device
		char buf[256];
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_NAME,
				sizeof(char[256]),
				&buf,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "OpenCL Device %d: Name = %s\n", i, buf);

        // and now we get the number of compute units of this device and printf it
		cl_uint units = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_MAX_COMPUTE_UNITS,
				sizeof(cl_uint),
				&units,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "OpenCL Device %d: Compute units = %u\n", i, units);

		// now we query and printf the largest work group of this device
		size_t gsize = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_MAX_WORK_GROUP_SIZE,
				sizeof(size_t),
				&gsize,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "OpenCL Device %d: Max. work group size = %d\n", i, (unsigned int)gsize);

		// finally, we query and printf the amount of constant memory (the fastest) the device has
		size_t maxConstantMemorySize = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
				sizeof(size_t),
				&maxConstantMemorySize,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "OpenCL Device %d: Max. const. mem. size = %d\n", i, (unsigned int)maxConstantMemorySize);	}

	if (!deviceFound) {
		fprintf(stderr, "Unable to select an appropriate device\n");
		exit(-1);
	}


	// creating the context
	// --------------------

	cl_context_properties cps[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties) platform,
		0
	};

	cl_context_properties *cprops = (NULL == platform) ? NULL : cps;
	// creates a context, which is a configured device and a queue of commands
	// properties (platform we're using etc.), num_devices (how many devices we want to use),
	// devices (array of devices to be used), pfn_notify (error callback), user_data, errcode_ret
	context = clCreateContext(
			cprops,
			1,
			&selectedDevice,
			NULL,
			NULL,
			&status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to open OpenCL context\n");
		exit(-1);
	}

    // checks if the context was properly created by querying the devices
	size_t deviceListSize;
    status = clGetContextInfo(
            context,
            CL_CONTEXT_DEVICES,
            32,
            devices,
            &deviceListSize);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to get OpenCL context info: %d\n", status);
		exit(-1);
    }

	// printfs the devices being used in the context we've just created
	for (i = 0; i < deviceListSize / sizeof(cl_device_id); ++i) {
		cl_device_type type = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_TYPE,
				sizeof(cl_device_type),
				&type,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		char *stype;
		switch (type) {
			case CL_DEVICE_TYPE_ALL:
				stype = "TYPE_ALL";
				break;
			case CL_DEVICE_TYPE_DEFAULT:
				stype = "TYPE_DEFAULT";
				break;
			case CL_DEVICE_TYPE_CPU:
				stype = "TYPE_CPU";
				break;
			case CL_DEVICE_TYPE_GPU:
				stype = "TYPE_GPU";
				break;
			default:
				stype = "TYPE_UNKNOWN";
				break;
		}
		fprintf(stderr, "[SELECTED] OpenCL Device %d: Type = %s\n", i, stype);

		char buf[256];
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_NAME,
				sizeof(char[256]),
				&buf,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "[SELECTED] OpenCL Device %d: Name = %s\n", i, buf);

		cl_uint units = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_MAX_COMPUTE_UNITS,
				sizeof(cl_uint),
				&units,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "[SELECTED] OpenCL Device %d: Compute units = %u\n", i, units);

		size_t gsize = 0;
		status = clGetDeviceInfo(devices[i],
				CL_DEVICE_MAX_WORK_GROUP_SIZE,
				sizeof(size_t),
				&gsize,
				NULL);
		if (status != CL_SUCCESS) {
			fprintf(stderr, "Failed to get OpenCL device info: %d\n", status);
			exit(-1);
		}

		fprintf(stderr, "[SELECTED] OpenCL Device %d: Max. work group size = %d\n", i, (unsigned int)gsize);
	}

	// creates a command queue and attaches it to our context
	// context, device (the device that will execute this queue), properties (), errcode_ret
	commandQueue = clCreateCommandQueue(
			context,
			devices[0],
			(cl_command_queue_properties)NULL,
			&status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL command queue: %d\n", status);
		exit(-1);
    }



    // creating buffers for our objects
    // --------------------------------

    // allocates space on the device memory to hold the objects that compose the scene
    // context (id), flags (rw, r, w etc), size (bytes of the buffer to be allocated),
    //      host_ptr (pointer to regular memory - a bit confusing...), errcode_ret
	objectBuffer = clCreateBuffer(
            context,
#ifdef __APPLE__
            CL_MEM_READ_WRITE, // NOTE: not READ_ONLY because of Apple's OpenCL bug
#else
			CL_MEM_READ_ONLY,
#endif
            sizeof(Object) * objectCount,
            NULL,
            &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL scene buffer: %d\n", status);
		exit(-1);
    }

    // now we issue a command on the queue that will write an array from regular memory (ram)
    // into the device's memory (vram)
    // command_queue (id), buffer (on which buffer to write), blocking_write (should block access
    //      to the src memory?), offset (from the src buffer), size (in bytes of src buffer),
    //      ptr (to the source buffer in ram), num_events_in_wait_list and event_wait_list
    //      (should we wait for these events to finish before exec'ing this command?),
    //      event (that represents this command - returned)
	status = clEnqueueWriteBuffer(
			commandQueue,
			objectBuffer,
			CL_TRUE,
			0,
			sizeof(Object) * objectCount,
			objects,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to write the OpenCL scene buffer: %d\n", status);
		exit(-1);
	}


    // creates a buffer to hold information regarding the camera (eye, target)
	cameraBuffer = clCreateBuffer(
            context,
#ifdef __APPLE__
            CL_MEM_READ_WRITE, // NOTE: not READ_ONLY because of Apple's OpenCL bug
#else
			CL_MEM_READ_ONLY,
#endif
            sizeof(Camera),
            NULL,
            &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL camera buffer: %d\n", status);
		exit(-1);
    }
	status = clEnqueueWriteBuffer(
			commandQueue,
			cameraBuffer,
			CL_TRUE,
			0,
			sizeof(Camera),
			&camera,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to write the OpenCL camera buffer: %d\n", status);
		exit(-1);
	}

	// allocates some buffers to hold the output from the opencl program
	allocateOutputBuffers();


	// creating the program from its source
	// ------------------------------------

	const char *sources = readFileByLine(kernelFileName);
	// context, count (number of lines in source), strings (array of lines),
	//      lengths (array of line lengths or NULL, if they're null-terminated), errcode_ret
	program = clCreateProgramWithSource(
        context,
        1,
        &sources,
        NULL,
        &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to open OpenCL kernel sources: %d\n", status);
		exit(-1);
    }

#ifdef __APPLE__
	status = clBuildProgram(program, 1, devices, "-I. -D__APPLE__", NULL, NULL);
#else
    // builds (compiles and links) the opencl program
    // program (id from createProgram), num_devices (devices that will run it),
    //      device_list, options ('command line'), ptn_notify (error cllbck), user_data
	status = clBuildProgram(program, 1, devices, "-I. ", NULL, NULL);
#endif
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to build OpenCL kernel: %d\n", status);

        size_t retValSize;
		status = clGetProgramBuildInfo(
				program,
				devices[0],
				CL_PROGRAM_BUILD_LOG,
				0,
				NULL,
				&retValSize);
        if (status != CL_SUCCESS) {
            fprintf(stderr, "Failed to get OpenCL kernel info size: %d\n", status);
			exit(-1);
		}

        char *buildLog = (char *)malloc(retValSize + 1);
        status = clGetProgramBuildInfo(
				program,
				devices[0],
				CL_PROGRAM_BUILD_LOG,
				retValSize,
				buildLog,
				NULL);
		if (status != CL_SUCCESS) {
            fprintf(stderr, "Failed to get OpenCL kernel info: %d\n", status);
			exit(-1);
		}
        buildLog[retValSize] = '\0';

		fprintf(stderr, "OpenCL Programm Build Log: %s\n", buildLog);
		exit(-1);
    }

    // creates our main kernel from the __kernel function called 'radianceGPU'
    // program (id), kernel_name (name of the __kernel function), errcode_ret
	kernel = clCreateKernel(program, "radianceGPU", &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to create OpenCL kernel: %d\n", status);
		exit(-1);
    }

	// determines the best work group size to use based on what the opencl driver thinks best
	// for the specific kernel we're going to execute (eg, how many registers it uses)
	size_t gsize = 0;
	status = clGetKernelWorkGroupInfo(kernel,
			devices[0],
			CL_KERNEL_WORK_GROUP_SIZE,
			sizeof(size_t),
			&gsize,
			NULL);

	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to get OpenCL kernel work group size info: %d\n", status);
		exit(-1);
	}

	workGroupSize = (unsigned int) gsize;
	fprintf(stderr, "OpenCL Device 0: kernel work group size = %d\n", workGroupSize);

	if (forceWorkSize > 0) {
		fprintf(stderr, "OpenCL Device 0: forced kernel work group size = %d\n", forceWorkSize);
		workGroupSize = forceWorkSize;
	}
}

///
/// Asks the device to execute our kernel, properly configured considering the work-groups and items.
///
static void executeKernel() {
    // globalThreads is the number of workGroups we need to compute all pixels 1 time
	size_t globalThreads[1];
	globalThreads[0] = width * height;
	if (globalThreads[0] % workGroupSize != 0) {
		globalThreads[0] = (globalThreads[0] / workGroupSize + 1) * workGroupSize;
        // suppose: 800x600 (ie, 480.000) and 1024 workGroupSize
        // numberOfWorkGroups = (479.700) * 1024 == 479.232	workgroups
	}
	// localThreads is the number of work-items a work-group has
	size_t localThreads[1];
	localThreads[0] = workGroupSize;

	// adds a command to execute our kernel
	// command_queue (id), kernel (id), work_dim (number of dimensions we'll provide the work-items),
	//      global_work_offset, global_work_size (number of global work-groups we need),
	//      local_work_size (how many work-items make a work-group - ), events...
	cl_int status = clEnqueueNDRangeKernel(
			commandQueue,
			kernel,
			1,
			NULL,
			globalThreads,
			localThreads,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to enqueue OpenCL work: %d\n", status);
		exit(-1);
	}
}

#define SET_KERNEL_ARG(kernel, i, type, value) do {\
        cl_int status = clSetKernelArg(kernel, i, sizeof(type), (void*)(value));\
        if (status != CL_SUCCESS) {\
            fprintf(stderr, "Failed to set OpenCL arg. #%d: %d\n", (i+1), status);\
        }\
    } while(0)

///
/// Asks the device to draw as many times as it can in 500ms by calling executeKernel() with the
/// updated values. Then, updates the color (output) buffer from vram to ram so glut redraws an
/// updated color buffer.
///
void updateRendering() {
	double frameStartTime = wallClockTime();
	int startSampleCount = currentSample;

    // sets the arguments for the kernel (the color buffer)
    // kernel (id), arg_index (index of this argument in the kernel function),
    //      arg_size (bytes of this argument), arg_value
    SET_KERNEL_ARG(kernel, 0, cl_mem, &colorBuffer);

	// now, sets the value for the random seed buffer
    SET_KERNEL_ARG(kernel, 1, cl_mem, &seedBuffer);

	// now, for the objects buffer
    SET_KERNEL_ARG(kernel, 2, cl_mem, &objectBuffer);

	// and the camera buffer
    SET_KERNEL_ARG(kernel, 3, cl_mem, &cameraBuffer);

	// and the objects count
    SET_KERNEL_ARG(kernel, 4, unsigned int, &objectCount);

	// and the scene width
    SET_KERNEL_ARG(kernel, 5, int, &width);

	// and height
    SET_KERNEL_ARG(kernel, 6, int, &height);

	// sets the value of the currentSample (what is this) in the kernel
    SET_KERNEL_ARG(kernel, 7, int, &currentSample);

	// sets the value of the pixel buffer (what is being drawn by opengl, gamma-corrected and in 4 bytes)
    SET_KERNEL_ARG(kernel, 8, cl_mem, &pixelBuffer);

	// sets the value of the debug buffer (to debug the Object struct alignment/packing)
    SET_KERNEL_ARG(kernel, 9, cl_mem, &debugBuffer);

    // sets the number of lights in the scene
    SET_KERNEL_ARG(kernel, 10, unsigned int, &lightCount);


	// asks the device to execute the kernel
	if (currentSample < 20) {
		executeKernel();
		currentSample++;
	} else {
		// if we are past the initial 20 samples, update the window less frequently
		// (every 500ms instead of for every sample)
		const float k = min(currentSample - 20, 100) / 100.f;
		const float tresholdTime = 0.5f * k;
		for (;;) {
			executeKernel();

			// waits for the queue to be fully executed (ie, become empty)
			clFinish(commandQueue);
			currentSample++;

			// keeps executing the kernel for more samples for at most 500ms, then proceed
			const float elapsedTime = wallClockTime() - frameStartTime;
			if (elapsedTime > tresholdTime)
				break;
		}
	}

	cl_int status;

	// now that we have executed the kernel as many times as we could on 500ms, we want to read
	// the results back from the opencl program (vram -> ram)
	// command_queue (id), buffer (memory address of what we want to read),
    //      blocking_read (should block the program here?), offset (on the src buffer),
    //      size (in bytes of the data being read), ptr (pointer to dest in ram), events...
	status = clEnqueueReadBuffer(
			commandQueue,
			pixelBuffer,
			CL_TRUE,
			0,
			width * height * sizeof(unsigned int),
			pixels,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to read the OpenCL pixel buffer: %d\n", status);
		exit(-1);
	}

	status = clEnqueueReadBuffer(
			commandQueue,
			debugBuffer,
			CL_TRUE,
			0,
			10 * sizeof(int),
			debug,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to read the OpenCL debug buffer: %d\n", status);
		exit(-1);
	} else {
	    // here we can show any debug info from the kernel... eg, those from the debugBuffer
        //printf("sizeof(Object) [CPU]: %d\n", sizeof(Object));
        //printf("sizeof(Object) [GPU]: %d\n", debug[0]);
	}


	updateRenderingStatistics(frameStartTime, startSampleCount);
}

/// Calculates the current rendering statistics: ellapsed time (beginning and this frame),
/// total passes, samples per second
///
/// Rendering time: total (this frame)
/// Passes: total (per second)
void updateRenderingStatistics(double frameStartTime, int startSampleCount) {
	// time spent on this "glut frame"... updates the string samples/sec to render it
	const double elapsedTimeThisFrame = wallClockTime() - frameStartTime;
	const double elapsedTime = wallClockTime() - startRenderingTime;
	char timeSinceBeginning[30];
	getHumanReadableTime(elapsedTime, timeSinceBeginning);
	sprintf(captionLine1, "Time:   %s  (%.3fs frames/s)", timeSinceBeginning, elapsedTimeThisFrame);

	const int samples = currentSample - startSampleCount;
	const double sampleSec = samples * height * width / elapsedTimeThisFrame;
	sprintf(captionLine2, "Passes: %d  (%.1fK samples/s)", currentSample, sampleSec / 1000.f);
}

///
/// Reinitializes the objects in the scene, called when the user changes them.
///
void reInitSceneObjects() {
	currentSample = 0;

	// writes the scene objects buffer to vram again
    cl_int status = clEnqueueWriteBuffer(
			commandQueue,
			objectBuffer,
			CL_TRUE,
			0,
			sizeof(Object) * objectCount,
			objects,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to write the OpenCL scene buffer: %d\n", status);
		exit(-1);
	}
}

///
/// Reinitializes rendering the scene camera, called when the viewpoint (camera) has changed.
///
void reInitViewpointDependentBuffers(const int reallocBuffers) {
    // when ' ' is pressed we want to throw away all our buffers and realloc them
	if (reallocBuffers) {
		freeBuffers();
		updateCameraBasis(&camera);
		allocateOutputBuffers();
	} else {
		updateCameraBasis(&camera);
	}

	// writes the camera buffer back to vram, as we've just changed it on the ram
	cl_int status = clEnqueueWriteBuffer(
			commandQueue,
			cameraBuffer,
			CL_TRUE,
			0,
			sizeof(Camera),
			&camera,
			0,
			NULL,
			NULL);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Failed to write the OpenCL camera buffer: %d\n", status);
		exit(-1);
	}

	currentSample = 0;
	startRenderingTime = wallClockTime();
}

///
/// Main function. Checks for cmd arguments, sets up the camera, opencl and glut.
/// Then delivers to glut (glutMainLoop).
///
int main(int argc, char *argv[]) {
	fprintf(stderr, "Usage: %s\n", argv[0]);
	fprintf(stderr, "Usage: %s <use CPU/GPU device (0=CPU or 1=GPU)> <workgroup size (0=default value or anything > 0 and power of 2)> <kernel file name> <window width> <window height> <scene file>\n", argv[0]);

	char* sceneName;
	if (argc == 7) {
		useGPU = atoi(argv[1]);
		forceWorkSize = atoi(argv[2]);
		kernelFileName = argv[3];
		width = atoi(argv[4]);
		height = atoi(argv[5]);
		sceneName = argv[6];
	} else if (argc == 1) {
		useGPU = 1;
		forceWorkSize = 0;
		kernelFileName = "path-tracing.cl";
		width = 480;
		height = 320;
		sceneName = "scenes/cornell.txt";
	} else {
		exit(-1);
    }

//    const int numberOfObjects = readScene(sceneName);
    readScene(sceneName);

	updateCameraBasis(&camera);

    setenv("CUDA_CACHE_DISABLE", "1", 1);
	setUpOpenCL();

	char windowTitle[150];
	sprintf(windowTitle, "FegemoPT: %s (%d objects)", sceneName, objectCount);
	initGlut(argc, argv, windowTitle);

    glutMainLoop();

	return 0;
}
