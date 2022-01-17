//used for .cu file
#include <typeinfo>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include "manager.h"
#include <cusparse_v2.h>
#include "cuda_runtime_api.h"

unordered_map<int,vector<vector<double>>> recordmapAA;
unordered_map<int,vector<vector<double>>> recordmapBB;

 

//kernel function
//__global__ void compute_1D_1D(int* d_out_ix, int* d_out_iy, double* d_in_ra1, double* d_in_dec1, double* d_in_ra2, double* d_in_dec2, double* d_out_dis, unsigned int nx, unsigned int ny){
//
//	unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
//	
//	if (ix+1 > nx) return;
//	
//	double d1 = d_in_ra1[ix], d2 = d_in_dec1[ix];
//	
//	for (unsigned int iy = 0; iy < ny; ++iy){
//		double tmp1 = d1 - d_in_ra2[iy];
//		double tmp2 = d2 - d_in_dec2[iy];
//		double tmp3 = d2 + d_in_dec2[iy];
//		tmp2 = tmp2*tmp2;
//		tmp3 = cos(tmp3/360.0*PI); // try to define 1/360 as macro
//		tmp1 = tmp1 * tmp3;
//		tmp1 = tmp1 * tmp1;
//		tmp1 = tmp1 + tmp2; // tmp1 is distance
//		//d_out_dis[(iy)*(nx)+ix] = tmp1;
//    if (tmp1<DIS){
//      size_t index = iy*nx+ix>>5; //相当于/32
//      size_t position = (iy*nx+ix) %32;
//      atomicOr(&d_out_ix[index], 1<<position-1);
//      atomicOr(&d_out_iy[index], 1<<position-1);
//    }
//    
//    //if (tmp1<DIS){
//      // d_out_dis[(iy)*(nx)+ix] = 1;
//   // }
//		
//	}
//	
//	return;
//}


__device__ double compute_distance(double rax,double decx,double ray, double decy){
    double tmp1 = rax - ray;
    double tmp2 = decx - decy;
    double tmp3 = decx + decy;
    tmp2 = tmp2*tmp2;
    tmp3 = cos(tmp3/360.0*PI); 
    tmp1 = tmp1 * tmp3;
    tmp1 = tmp1 * tmp1;
    tmp1 = tmp1 + tmp2; 
    return sqrt(tmp1);
}

__device__ int AtomicAdd(unsigned int *address, int incr)
{
    // Create an initial guess for the value stored at *address.
    unsigned int guess = *address;
    unsigned int oldValue = atomicCAS(address, guess, guess + incr);

    // Loop while the guess is incorrect.
    while (oldValue != guess)
    {
        guess = oldValue;
        oldValue = atomicCAS(address, guess, guess + incr);
    }

    return oldValue;
}

__global__ void compute_1D_1D(double* d_in_ra1, double* d_in_dec1, double* d_in_ra2, double* d_in_dec2, unsigned int nx, unsigned int ny, int  *ixl, int *iyl,unsigned int *sharedInteger){
   
   //__shared__ unsigned int *sharedInteger;
	unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
	if (ix+1 > nx) return;
  	double d1 = d_in_ra1[ix], d2 = d_in_dec1[ix];
  	
	for (unsigned int iy = 0; iy < ny; ++iy){
     double distance = compute_distance(d1,d2,d_in_ra2[iy],d_in_dec2[iy]);
  	
      if (distance < DIS){
        
         //atomicAdd(&(temp[buffer[i]]),1);
//        //d_out_dis[(iy)*(nx)+ix] = distance;
//       // printf("-----::%d\n",*sharedInteger);
        atomicAdd(sharedInteger, 1); 
//        ixl[*sharedInteger] = ix;
//        iyl[*sharedInteger] = iy; 
        //__syncthreads();

      }
  	}
  
	return;
}



//typedef struct seq_s seq_t;
//struct seq_s {
//    unsigned int* ixlist;
//    unsigned int* iylist;
//};


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


int main(int argc, char **argv){
    vector<double> gpu_info = gpu_helper();
    double total_db = gpu_info[0];
    double free_db = gpu_info[1];
 
     
    double AllStart = cpuSecond();
    double iStart,iElaps;
    //[CPU 1]: read the csv input into recordmap
    recordmapAA = read_to_unordered("data/gaiatest.csv",recordmapAA);
    iElaps = cpuSecond() - AllStart;
	  printf("[Time:] A read to unorded time is %f s\n", iElaps);
    recordmapBB = read_to_unordered("data/gaiatest2.csv",recordmapBB);
    iElaps = cpuSecond() - AllStart;
    printf("[Time:] A+B read to unorded time is %f s\n", iElaps);
    
    //[CPU 2]: get basic info of input reocrd A and B
    unordered_map<int,vector<vector<double>>>::iterator iterA;
    unordered_map<int,vector<vector<double>>>::iterator iterB;
	  int linenumA,linenumB;
  	for (iterA=recordmapAA.begin();iterA!=recordmapAA.end();iterA++){
	    //cout <<iterA->first << " " << iterA->second.size() <<endl;
	    linenumA=linenumA + iterA->second.size();
	  }
   	for (iterB=recordmapBB.begin();iterB!=recordmapBB.end();iterB++){
	    //cout <<iterB->first << " " << iterB->second.size() <<endl;
	    linenumB=linenumB + iterB->second.size();
	  }
     
    int averageA = linenumA/recordmapAA.size();
    int averageB = linenumB/recordmapBB.size();
    printf("[A info]:average number of lines in blocks : %d\n",averageA);
  	printf("[A info]:total number of lines in csv file: %d\n",linenumA);
  	printf("[A info]:total number of divided blocks: %d\n",recordmapAA.size());
  	printf("[B info]:average number of lines in blocks : %d\n",averageB);
  	printf("[B info]:total number of lines in csv file: %d\n",linenumB);
  	printf("[B info]:total number of divided blocks: %d\n",recordmapBB.size());
   
    //get index list and shared index list
    //vector<vector<int>> idlistA,idlistB;
    //idlistA =  get_idlist(recordmapAA,idlistA);
    //idlistB =  get_idlist(recordmapBB,idlistB);
    vector<int> sharedlist = get_shared_id(recordmapAA,recordmapBB);

    //[CPU 3*]: combine blocks
    //iStart = cpuSecond();
    //vector<vector<int>> combine;
    //vector <int> temp;
    //int threshold= (averageA+averageB)/2;
    //combine=combineblocks(threshold,combine,sharedlist,idlistA,idlistB,0,0,0,temp,0);
    //printf("[Debug:] the num of blocks after combined is: %d\n",combine.size());
    //iElaps = cpuSecond() - iStart;
    //printf("[Time:] time for combining is %f s\n", iElaps);
    
     
    //----------------------------------------------------------------------------------------
    //[GPU 1]: Initialization
//   	int dev = 0;
//	  CHECK(cudaSetDevice(dev));
     
     
    	// 3.2 initialize multi-streams
	int stream_size = 3;
//	printf("[DEBUG]The size of shared list [num of Stream] is %d\n",stream_size);
	cudaStream_t stream[stream_size];
	for(int i=0;i<stream_size;i++){
		cudaStreamCreateWithFlags(&stream[i],cudaStreamNonBlocking);
	}	
 
   // define variety for communication of CPU and GPU
    double *h_in_ra1, *h_in_dec1, *h_in_ra2, *h_in_dec2;
	  //double  *h_out_dis;
	  double *d_in_ra1, *d_in_dec1, *d_in_ra2, *d_in_dec2;
	  //double  *d_out_dis;
    int *d_out_ixl;
    int *d_out_iyl;
    int *h_out_ixl;
    int *h_out_iyl;
     
	  double *h_in_ra_x, *h_in_dec_x, *h_in_ra_y, *h_in_dec_y;
    //size_t nBytes, nBytes2, nBytes3;
 

  vector<vector<double>> matchresult;
//     int m=5;
//    int n=10000000;
//    vector<vector<double>> a(n);
//    //a.resize(n*sizeof(double));
//   for (int i=0;i<n;i++){
//      a[i].resize(m);
//   }
   
  // 3.4 malloc device memory space paged locked memory
	
	CHECK(cudaMallocHost((double**)&h_in_ra1, N*sizeof(double)));
	CHECK(cudaMallocHost((double**)&h_in_ra2, N*sizeof(double)));
	CHECK(cudaMallocHost((double**)&h_in_dec1, N*sizeof(double)));
	CHECK(cudaMallocHost((double**)&h_in_dec2, N*sizeof(double)));
	//CHECK(cudaMallocHost((double**)&h_out_dis, BLOCK_MAX_X*BLOCK_MAX_Y*sizeof(double)));
	//CHECK(cudaMallocHost((double**)&h_out_dis,  BLOCK_MAX_X*BLOCK_MAX_Y*sizeof(double)));
  //CHECK(cudaMallocHost((int**)&h_out_ix_iy,  BLOCK_MAX_X*sizeof(int*)));
  CHECK(cudaMallocHost((int**)&h_out_ixl,  BLOCK_MAX_X*sizeof(int*)));
  CHECK(cudaMallocHost((int**)&h_out_iyl,  BLOCK_MAX_Y*sizeof(int*)));
  
   //h_out_ix_iy = (int**)malloc(BLOCK_MAX_X*BLOCK_MAX_Y*sizeof(int*));
 
  
	// 4 scanf data from filesystem
	// loop arguments
  unsigned int resultcount = 0;
	unsigned int nx, ny;
	int s_index=0;

//得到 inputdata and nxylist
    iStart = cpuSecond();
    vector<vector<double*>> inputdata;
    vector<vector<int>> nxylist;
    for(int i=0; i<sharedlist.size(); i++){
      vector <double*> subinput;
      vector <int> subnxy;
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
      subinput.push_back(h_in_ra_x);
      subinput.push_back(h_in_dec_x);
      subinput.push_back(h_in_ra_y);
      subinput.push_back(h_in_dec_y);
      inputdata.push_back(subinput);
      subinput.clear();
      
      subnxy.push_back(nx);
      subnxy.push_back(ny);
      nxylist.push_back(subnxy);
      subnxy.clear();

  }
  iElaps = cpuSecond() - iStart;
  printf("[Info] input data generation cost %f s\n", iElaps);
    
   int small_bloack_count=0;
  // cuda memory malloc
  //得到 offset band
  vector<vector<unsigned int>> offset_band;
  for(unsigned int i=0; i<sharedlist.size(); i++){
    nx = nxylist[i][0];
    ny = nxylist[i][1];
    //对每一个二次划分之后的小计算块
		for (unsigned int data_x_offset = 0; data_x_offset < nx; data_x_offset += BLOCK_MAX_X){
			for (unsigned int data_y_offset = 0; data_y_offset < ny; data_y_offset += BLOCK_MAX_Y){
        vector<unsigned int> tempband;
				unsigned int data_x_band = min(nx-data_x_offset, BLOCK_MAX_X), 
				data_y_band = min(ny-data_y_offset, BLOCK_MAX_Y);
        tempband.push_back(i);
			  tempband.push_back(data_x_band);
        tempband.push_back(data_y_band);
        tempband.push_back(data_x_offset);
        tempband.push_back(data_y_offset);
        offset_band.push_back(tempband);
        tempband.clear();
        small_bloack_count+=1;
       }
     }
   }
   
   printf("[Debug:]the total small bloack num is: %d\n",small_bloack_count);
   printf("[Debug:]the total small bloack num is: %d\n",offset_band.size());
 
   //分配cuda大小
   iStart = cpuSecond();
   vector<size_t> maxnbytes;
   maxnbytes = get_cudamalloc_size(offset_band);
 	 CHECK(cudaMalloc((void **) &d_in_ra1, maxnbytes[0]));
   CHECK(cudaMalloc((void **) &d_in_dec1, maxnbytes[0]));
   CHECK(cudaMalloc((void **) &d_in_ra2, maxnbytes[1]));
   CHECK(cudaMalloc((void **) &d_in_dec2,maxnbytes[1]));
   //CHECK(cudaMalloc((void **) &d_out_dis, maxnbytes[2]));
   //CHECK(cudaMalloc((void **) &d_out_ix_iy,maxnbytes[3]));
   CHECK(cudaMalloc((void **) &d_out_ixl, maxnbytes[3]));
   CHECK(cudaMalloc((void **) &d_out_iyl, maxnbytes[4]));
   
   
   iElaps = cpuSecond() - iStart;
   printf("[Info] cuda Malloc cost %f s\n", iElaps);
   
   
   
   //!!!!!
//   seq_s *seqs;
//	 seqs = (seq_s *)malloc(sizeof(seq_s));	 
//	 seqs->ixlist = (unsigned int *)malloc(sizeof(int));
//   seqs->iylist = (unsigned int *)malloc(sizeof(int));
//   
// 	 seq_s * tmp_seq = (seq_s *)malloc(sizeof(seq_s));
//   memcpy(&tmp_seq, seqs,  sizeof(seq_s)); 
//  
//   cudaMalloc((void **) &(tmp_seq->ixlist),  sizeof(int)); 
//   cudaMalloc((void **) &(tmp_seq->iylist), sizeof(int));
//   cudaMemcpy(tmp_seq->ixlist, seqs->ixlist, sizeof(int), cudaMemcpyHostToDevice);
//   cudaMemcpy(tmp_seq->iylist, seqs->iylist, sizeof(int), cudaMemcpyHostToDevice);
//   
//   
// 	 seq_s * gpu_seq;
//	 cudaMalloc((void**)&gpu_seq,sizeof(seq_s));
//	 cudaMemcpy(gpu_seq, tmp_seq, sizeof(seq_s), cudaMemcpyHostToDevice);
    //!!!!!
   
   
   

   int checknum=0;
   //calculate the distance  相当于是对每一个能放进GPU 全局内存的小计算块 1011个一共
   for (int m=0; m<offset_band.size(); m++){
        unsigned int *h_sharedInteger;
        CHECK(cudaMallocHost((int**)&h_sharedInteger, sizeof(int)));
        unsigned int *d_sharedInteger;
        CHECK(cudaMalloc((void **)&d_sharedInteger, sizeof(int)));
        CHECK(cudaMemset(d_sharedInteger, 0, sizeof(int)));
   
        int fileindex = offset_band[m][0];
        unsigned int data_x_band = offset_band[m][1];
        unsigned int data_y_band = offset_band[m][2];
        unsigned int data_x_offset = offset_band[m][3];
        unsigned int data_y_offset = offset_band[m][4];     
        int dimx = 256;
				int dimy = 1;
        if (argc > 1) dimx = atoi(argv[1]);
				dim3 block(dimx, dimy);
				dim3 grid((data_x_band + block.x - 1) / block.x, 1);
        
        // 6.2 execute the kernel and receive the result of GPU
				iStart = cpuSecond();
        int j = m % 3;                     
				CHECK(cudaMemcpyAsync(d_in_ra1, inputdata[fileindex][0]+data_x_offset, maxnbytes[0], cudaMemcpyHostToDevice,stream[j]));
				CHECK(cudaMemcpyAsync(d_in_dec1, inputdata[fileindex][1]+data_x_offset, maxnbytes[0], cudaMemcpyHostToDevice,stream[j]));
				CHECK(cudaMemcpyAsync(d_in_ra2, inputdata[fileindex][2]+data_y_offset, maxnbytes[1], cudaMemcpyHostToDevice,stream[j]));
				CHECK(cudaMemcpyAsync(d_in_dec2, inputdata[fileindex][3]+data_y_offset,maxnbytes[1], cudaMemcpyHostToDevice,stream[j]));	
				//compute
        compute_1D_1D <<<grid, block,0,stream[j]>>>(d_in_ra1, d_in_dec1, d_in_ra2, d_in_dec2, data_x_band, data_y_band,d_out_ixl,d_out_iyl,d_sharedInteger);
        cudaStreamSynchronize(stream[j]); 
        iElaps = cpuSecond() - iStart;
				//printf("[Info] compute <<<(%d,%d), (%d,%d)>>> elapsed %f s\n", grid.x, grid.y, block.x, block.y, iElaps);
        CHECK(cudaMemcpyAsync(h_sharedInteger, d_sharedInteger, sizeof(int),cudaMemcpyDeviceToHost,stream[j]));
        // 
         //if (h_sharedInteger!=checknum){
           //printf("!!!!!!!%d\n", *h_sharedInteger);
           cout<<*h_sharedInteger<<endl;
           resultcount= resultcount + *h_sharedInteger;
        // }
        
         
         //CHECK(cudaMemcpyAsync(h_out_ixl, d_out_ixl,maxnbytes[3], cudaMemcpyDeviceToHost,stream[j]));
         //CHECK(cudaMemcpyAsync(h_out_iyl, d_out_iyl,maxnbytes[4], cudaMemcpyDeviceToHost,stream[j]));
   
//        //CHECK(cudaMemcpyAsync(h_out_dis, d_out_dis, maxnbytes[2], cudaMemcpyDeviceToHost,stream[j]));
//        CHECK(cudaMemcpyAsync(h_out_dis, d_out_dis, maxnbytes[2], cudaMemcpyDeviceToHost,stream[j]));
       
       //checknum=h_sharedInteger;
       
        CHECK(cudaFree(d_sharedInteger));
        CHECK(cudaFreeHost(h_sharedInteger));
   }
   
   
//   for(int i=0;i<3;i++){
//     cudaStreamSynchronize(stream[i]); 
//   }
//    
         
    //printf("gqqq%ld\n",maxnbytes[3]);
   //CHECK(cudaMemcpy(&h_sharedInteger, d_sharedInteger, sizeof(int),cudaMemcpyDeviceToHost));
  //printf("4 x 128 increments led to value of %d\n", h_sharedInteger);
  printf("!!!!!!!--%ld--\n",resultcount);
   
   
  
   
//   int countcount=0;
//   // host out data --- output the crossmatch result  
//   //#pragma omp parallel for shared(a) private(m)     
//   for (int m=0; m<offset_band.size(); m++){
//      //printf("[DEBUG0] lalala %d\n: ",countcount);
//      unsigned int data_x_band = offset_band[m][1];
//      unsigned int data_y_band = offset_band[m][2];
//      unsigned int data_x_offset = offset_band[m][3];
//      unsigned int data_y_offset = offset_band[m][4];
//      iStart = cpuSecond();
//      //vector<double> tempresult;
//      	for (unsigned int ty = 0; ty<data_y_band; ty++){
//					for (unsigned int tx = 0; tx<data_x_band; tx++){
//						unsigned int tidx = (ty)*(data_x_band)+(tx);
//            //printf("[DEBUG] lalala %d\n: ",countcount);
//						if (h_out_dis[tidx] ==1){
//                                
//              a[countcount][0]=h_in_ra_x[tx+data_x_offset];
//              a[countcount][1]=h_in_dec_x[tx+data_x_offset];
//              a[countcount][2]=h_in_ra_y[ty+data_y_offset];
//              a[countcount][3]=h_in_dec_y[ty+data_y_offset];
//              a[countcount][4]=sqrt(h_out_dis[tidx]);
//              countcount=countcount+1;
//              
//						}
//             //printf("[DEBUG] lalala %d\n: ",countcount);
//					}
//            //printf("[DEBUG] lalala %d\n: ",countcount);
//				}
//     	
//    	iElaps = cpuSecond() - iStart;
//    	printf("[Info] result generation cost %f s\n", iElaps);
//     //printf("[DEBUG] yoyoyo %d\n: ",m);
//     
//    }
    
    
    
    
   	//8 destroy memory
  	CHECK(cudaFreeHost(h_in_dec1));
  	CHECK(cudaFreeHost(h_in_dec2));
  	CHECK(cudaFreeHost(h_in_ra1));
  	CHECK(cudaFreeHost(h_in_ra2));
  	
   	CHECK(cudaFree(d_in_ra1));
  	CHECK(cudaFree(d_in_ra2));
  	CHECK(cudaFree(d_in_dec1));
  	CHECK(cudaFree(d_in_dec2));
   
   // 7 destroy streams
  	for(int i=0;i<stream_size;i++){
  		cudaStreamDestroy(stream[i]);
  	}
 
//  cout<<a[4000000][0]<<endl;
//  cout<<a[2000000][3]<<endl;
//  cout<<a[4500000][2]<<endl;
//  printf("[!!!] %d \n",countcount);
//  printf("[!!!] %.16lf \n",a[4000000][0]);
//  printf("[!!!] the match result size is %.16lf\n",a[2000000][3]);
//  printf("[!!!] the match result size is %.16lf\n",a[4500000][2]);
  
  
  //printf("[!!!] the match result size is %d\n",a.size());
  
  
  
  //CHECK(cudaFreeHost(h_out_dis));
	//CHECK(cudaFree(d_out_dis));
  
  CHECK(cudaFreeHost(h_out_ixl));
  CHECK(cudaFreeHost(h_out_iyl));
  CHECK(cudaFree(d_out_ixl));
  CHECK(cudaFree(d_out_iyl));

	//print all time
	iElaps = cpuSecond() - AllStart;
	printf("[Info]All time is %f s\n", iElaps);	   
   

return 0;

}