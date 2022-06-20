#ifndef HELPER_H
#define HELPER_H
#include "constant.h"


#define CHECK(call)						\
{								\
    const cudaError_t error = call;				\
    if (error != cudaSuccess){					\
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);	\
	fprintf(stderr, "code: %d, reason: %s\n", error,	\
		cudaGetErrorString(error));			\
    }								\
}



double cpuSecond();

float string_to_float(string str);

void file_to_string(vector<string> &record, const string& line, char delimiter);

double distance(double rax,double decx,double ray, double decy);


#endif  // HELPER_H