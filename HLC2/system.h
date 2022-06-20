#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "cuda_runtime_api.h"
#include <vector>
using namespace std;

double cpu_helper();

vector<double> gpu_helper();