#include <typeinfo>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include "manager.h"
#include <cusparse_v2.h>

unordered_map<int,vector<vector<double>>> recordmapAA;
unordered_map<int,vector<vector<double>>> recordmapBB;

__device__ double compute_distance(double rax,double decx,double ray, double decy){
    double tmp1 = rax - ray;
    double tmp2 = decx - decy;
    double tmp3 = decx + decy;
    tmp2 = tmp2*tmp2;
    tmp3 = cos(tmp3/360.0*PI); // try to define 1/360 as macro
    tmp1 = tmp1 * tmp3;
    tmp1 = tmp1 * tmp1;
    tmp1 = tmp1 + tmp2; 
    return tmp1;
}


//kernel function
__global__ void compute_1D_1D(double* d_in_ra1, double* d_in_dec1, double* d_in_ra2, double* d_in_dec2, int* d_out_dis, unsigned int nx, unsigned int ny,unsigned int *resultcount){

    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;	
    if (ix+1 > nx) return;	
	double d1 = d_in_ra1[ix], d2 = d_in_dec1[ix];
	for (unsigned int iy = 0; iy < ny; ++iy){
	    double tmp1 = compute_distance(d1,d2,d_in_ra2[iy], d_in_dec2[iy]);
	    if (tmp1<DIS){
	     d_out_dis[(iy)*(nx)+ix] = 1;
	     atomicAdd(resultcount, 1);
             }
	}
  
  //printf("max index: %d\n",maxindex);
	return;
}


int main(int argc, char **argv){

    double AllStart = cpuSecond();
    double iStart, iElaps;
  
    vector<double> parameters = parameter_decided();
    const int N = parameters[0];
    const int BLOCK_MAX_X = parameters[1];
    const int BLOCK_MAX_Y = parameters[2];
  
    //read the csv input into recordmap
    recordmapAA = read_to_unordered("data/twomasstest5.csv",recordmapAA);
    iElaps = cpuSecond() - AllStart;
    printf("A read to unorded time is %f s\n", iElaps);
    recordmapBB = read_to_unordered("data/twomasstest5.csv",recordmapBB);
    iElaps = cpuSecond() - AllStart;
    printf("A+B read to unorded time is %f s\n", iElaps);
    unordered_map<int,vector<vector<double>>>::iterator iterA;
    unordered_map<int,vector<vector<double>>>::iterator iterB;
    int linenumA,linenumB;
  	
   
    //get shared index list of recordmap A and B
    vector<int> sharedlist;
    sharedlist = get_shared_id(recordmapAA,recordmapBB);
    printf("shared index num: %d \n",sharedlist.size());
    iElaps = cpuSecond() - AllStart;
    printf("read csv file time is %f s\n", iElaps);
    //printf ("__________\n");
 
    vector<vector<double>> matchresult;
   
    // GPU Initialization
    // set up device
    int dev = 0;
    CHECK(cudaSetDevice(dev));
   
   
   // define variety for communication of CPU and GPU
     double *h_in_ra1, *h_in_dec1, *h_in_ra2, *h_in_dec2;
     double *h_out_ra1, *h_out_dec1, *h_out_ra2, *h_out_dec2;
   
     double *d_in_ra1, *d_in_dec1, *d_in_ra2, *d_in_dec2;
     double *d_out_ra1, *d_out_dec1, *d_out_ra2, *d_out_dec2;
    
     double *h_in_ra_x, *h_in_dec_x, *h_in_ra_y, *h_in_dec_y;
     int *h_out_dis,*d_out_dis;
     unsigned int *h_count, *d_count;

     CHECK(cudaMallocHost((int**)&h_count, sizeof(int)));
     CHECK(cudaMalloc((void **)&d_count, sizeof(int)));
     CHECK(cudaMemset(d_count, 0, sizeof(int)));
  
     // malloc device memory space
     CHECK(cudaMallocHost((double**)&h_in_ra1, N*sizeof(double)));
     CHECK(cudaMallocHost((double**)&h_in_ra2, N*sizeof(double)));
     CHECK(cudaMallocHost((double**)&h_in_dec1, N*sizeof(double)));
     CHECK(cudaMallocHost((double**)&h_in_dec2, N*sizeof(double)));
   
     CHECK(cudaMallocHost((int**)&h_out_dis, BLOCK_MAX_X*BLOCK_MAX_Y*sizeof(int)));
     CHECK(cudaMalloc((void **) &d_out_dis, BLOCK_MAX_X*BLOCK_MAX_Y*sizeof(int)));
    
    size_t nBytes, nBytes2, nBytes3;
    unsigned int resultcount = 0;
    unsigned int compress_resultcount = 0;
    unsigned int nx, ny;
    
    int streamcount=0;

     for(int i=0; i<sharedlist.size(); i++){
        vector<vector<double>> valuesA = recordmapAA[sharedlist[i]];
        vector<vector<double>> valuesB = recordmapBB[sharedlist[i]];
        for(int r=0; r<valuesA.size(); r++){
          h_in_ra1[r] = valuesA[r][2];
          h_in_dec1[r] = valuesA[r][3];
        }
        
        for(int t=0; t<valuesB.size(); t++){
          h_in_ra2[t] = valuesB[t][2];
          h_in_dec2[t] = valuesB[t][3];
        }
      
        // make it convenient to use varieties
    	if (valuesA.size() > valuesB.size()){
    		h_in_ra_x = h_in_ra1;
    		h_in_dec_x = h_in_dec1;
    		h_in_ra_y = h_in_ra2;
    		h_in_dec_y = h_in_dec2;
    		nx = valuesA.size();
    		ny = valuesB.size();
       	 }else{		
		h_in_ra_x = h_in_ra2;
    		h_in_dec_x = h_in_dec2;
    		h_in_ra_y = h_in_ra1;
    		h_in_dec_y = h_in_dec1;
    		nx = valuesB.size();
    		ny = valuesA.size();
    	}
       
      
      // for each calculation block
    	for (unsigned int data_x_offset = 0; data_x_offset < nx; data_x_offset += BLOCK_MAX_X){
    		for (unsigned int data_y_offset = 0; data_y_offset < ny; data_y_offset += BLOCK_MAX_Y){
            		streamcount = streamcount + 1; 
            
    			unsigned int data_x_band = min(nx-data_x_offset, BLOCK_MAX_X), 
    			data_y_band = min(ny-data_y_offset, BLOCK_MAX_Y);
    				
			unsigned int *h_sharedInteger,*d_sharedInteger;
			CHECK(cudaMallocHost((int**)&h_sharedInteger, sizeof(int)));
			CHECK(cudaMalloc((void **)&d_sharedInteger, sizeof(int)));
			CHECK(cudaMemset(d_sharedInteger, 0, sizeof(int)));
            
    			// dynamically set up the size of data
    			nBytes = (data_x_band)*sizeof(double);
    			nBytes2 = (data_y_band)*sizeof(double);
    			nBytes3 = (data_x_band)*(data_y_band)*sizeof(int);
    				
    			int dimx = 128;
    			int dimy = 1;
    
    			if (argc > 1) dimx = atoi(argv[1]);
    			dim3 block(dimx, dimy);
    			dim3 grid((data_x_band + block.x - 1) / block.x, 1);
    				
    			// execute the kernel and receive the result of GPU
    			iStart = cpuSecond();
    			CHECK(cudaMalloc((void **) &d_in_ra1, nBytes));
    			CHECK(cudaMalloc((void **) &d_in_dec1, nBytes));
    			CHECK(cudaMalloc((void **) &d_in_ra2, nBytes2));
    			CHECK(cudaMalloc((void **) &d_in_dec2,nBytes2));
    			//CHECK(cudaMalloc((void **) &d_out_dis, nBytes3));
    
          		 //int j = streamcount % 2;                                  
    			CHECK(cudaMemcpyAsync(d_in_ra1, h_in_ra_x+data_x_offset, nBytes, cudaMemcpyHostToDevice));
    			CHECK(cudaMemcpyAsync(d_in_dec1, h_in_dec_x+data_x_offset, nBytes, cudaMemcpyHostToDevice));
    			CHECK(cudaMemcpyAsync(d_in_ra2, h_in_ra_y+data_y_offset, nBytes2, cudaMemcpyHostToDevice));
    			CHECK(cudaMemcpyAsync(d_in_dec2, h_in_dec_y+data_y_offset,nBytes2, cudaMemcpyHostToDevice));	
             		
            		compute_1D_1D <<<grid, block>>>(d_in_ra1, d_in_dec1, d_in_ra2, d_in_dec2, d_out_dis, data_x_band, data_y_band,d_sharedInteger);
            
    			CHECK(cudaMemcpyAsync(h_out_dis, d_out_dis, nBytes3, cudaMemcpyDeviceToHost));
            
    			 CHECK(cudaMemcpy(h_sharedInteger, d_sharedInteger, sizeof(int),cudaMemcpyDeviceToHost));
 
			//cout<<*h_sharedInteger<<endl;         
			resultcount= resultcount + *h_sharedInteger; 
  
    			CHECK(cudaFree(d_in_ra1));
    			CHECK(cudaFree(d_in_ra2));
    			CHECK(cudaFree(d_in_dec1));
    			CHECK(cudaFree(d_in_dec2));
			CHECK(cudaFree(d_sharedInteger));
			CHECK(cudaFreeHost(h_sharedInteger));

    			}
    		}
       
       
    		iElaps = cpuSecond() - iStart;
    		//printf("[Info]File %d is done, elapsed %f s.\n", sharedlist[i], iElaps);
  	}
     

	 //check the correctness
	 printf("!!!!!!!--%ld--\n",resultcount);  
	 printf("[result:] %d \n",matchresult.size());
   
    
	//destroy memory
  	CHECK(cudaFreeHost(h_in_dec1));
  	CHECK(cudaFreeHost(h_in_dec2));
  	CHECK(cudaFreeHost(h_in_ra1));
  	CHECK(cudaFreeHost(h_in_ra2));
  	CHECK(cudaFreeHost(h_out_dis));
	CHECK(cudaFree(d_out_dis)); 
	CHECK(cudaFree(d_count));
	CHECK(cudaFreeHost(h_count));  
    

  	iElaps = cpuSecond() - AllStart;
  	printf("[Info]All time is %f s\n", iElaps);	
   
  	return 0;
  }
  
