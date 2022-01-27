#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include <set>

#include <sys/time.h>
#include "chealpix.h"
#include <queue>
#include <pthread.h>
#include <unordered_map>
#include <map>
#include <string>
#include <atomic> 


using namespace std;

#define R_A 0.001388889/5
#define R_B 0.001388889/5
#define RADIANS 57.29577951308232311024
#define D2R  (1.7453292519943295769e-2)
#define PI 3.1415926535
#define DIS 9.0*(R_A*R_A+R_B*R_B)

#define N (unsigned long)(1<<19) // 524,288
//#define SHARED_Y (1<<11) //2048
#define BLOCK_MAX_X (1<<15) //32,768
#define BLOCK_MAX_Y (1<<14) //16,384

