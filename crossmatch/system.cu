#include "system.h"
//source reading module
double cpu_helper()
{
    FILE *fp = NULL;
    FILE *fp2 = NULL;
    fp=popen("cat /proc/meminfo | grep MemTotal:|sed -e 's/.*:[^0-9]//'","r");
    fp2=popen("cat /proc/meminfo | grep MemAvailable:|sed -e 's/.*:[^0-9]//'","r");

    char buf[100];
    char buf2[100];
    double avaliableCPU;
    double totalCPU;
    while(memset(buf, 0, sizeof(buf)), fgets(buf, sizeof(buf) - 1, fp) != 0 ) {
      totalCPU = strtod(buf,NULL);
      //printf("%s", buf);
      //printf( "%f\n", totalCPU );
    }
    
    while(memset(buf2, 0, sizeof(buf)), fgets(buf2, sizeof(buf2) - 1, fp2) != 0 ) {
      //printf("%s", buf2);
      avaliableCPU = strtod(buf2,NULL);
      //printf( "%f\n", avaliableCPU);
    }

    pclose(fp);
    
    return avaliableCPU;
}

vector<double> gpu_helper()
{
	vector<double> gpu_info;
	size_t free_byte;
	size_t total_byte;
	cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);
	if (cudaSuccess != cuda_status) {
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
		exit(1);
	}

	double free_db = (double)free_byte;
	double total_db = (double)total_byte;
	//double used_db_1 = (total_db - free_db) / 1024.0 / 1024.0;
  	double total_db_1 = total_db /1024 /1024;
  	double free_db_1 = free_db /1024/1024;
  
  	//gpu_info.push_back(used_db_1);
  	gpu_info.push_back(total_db_1);
  	gpu_info.push_back(free_db_1);
  
  	return gpu_info;
}
