#include "manager.h"
#include "ang2pix_nest.c"
#include "mk_xy2pix.c"
#include "system.h"
#include "HTM.h"
#include "HTM.c"

int power (double gpumem){
  int coefficient=0; 
  #pragma omp parallel for
    for(int i=0;i<21;i++){
       double left = 1<<i;
       double right = 1<<(i+1);
       if (gpumem>left && gpumem <right) {
         coefficient = i+1;
       }
    }
  
  //cout<<"testtest:"<<coefficient<<endl;
  return coefficient;
}

vector<double> parameter_decided(){
    vector<double> parameter;
    vector<double> gpu_info = gpu_helper();
    double free_db,availableCPU;
    free_db = gpu_info[1];
   
    availableCPU = cpu_helper();
    int nnn;
    int left, right;
    //---------------cpu--------------------
    right = power(availableCPU/1024)+20;
    
    //nnn = (availableCPU * 1024) / 32 ; 

    int j=1;
    
    while(j<30){
        left =15+j; 
        if (left < right){
          j=j+1;
        }
        else{
          break;
      }
    
    }
    cout<<j<<endl;
    nnn= (1<<(j-1));
    parameter.push_back(nnn);
    
    
    //----------GPU---------------
    //parameter.push_back(free_db);
    
    int compare_left, compare_right;
    int compare_left2;
    int block_max_x, block_max_y;

    compare_right = power(free_db)+20;
    int i=1;
    while(i<20){
      compare_left = 8*4*(1<<i)+8*(1<<i)*(1<<i);
      compare_left2 = 3+i+i;
      if (compare_left2 < compare_right){
          i=i+1;
      }
      else{
          break;
      }
      
    }
    
    
    cout<<i<<endl;
    block_max_x= (1<<(i-1));
    block_max_y = block_max_x /2 ;
    
    parameter.push_back(block_max_x);
    parameter.push_back(block_max_y);
    
  return parameter;
}

//calculating Healpix index
int hp_nestid(double ra, double dec, int order)
{
	long id;
	long nside = pow(2, order);
	double theta = (90.0 - dec) * D2R;
	double phi = ra * D2R;
	ang2pix_nest(nside, theta, phi, &id);
	return id;
}


// read the star records from csv file to [format: unordered_map<hp index, corresponding star records> ]
unordered_map<int,vector<vector<double>>> read_to_unordered(string path,unordered_map<int,vector<vector<double>>> recordmap2){
	  vector<double>b;
    vector<string> row;
    string line;
    string filename;
    ifstream in(path);
    if (in.fail())  { cout << "File not found" <<endl; }
    
    while(getline(in, line)  && in.good() )
    {
        file_to_string(row, line, ',');  
        #pragma omp parallel for
        for(int i=0, leng=row.size(); i<leng; i++){
              b.push_back(string_to_float(row[i]));
        }
     
        int hpid;
        int htmid;
    	int q_num,q_num1,q_num2,q_num3,q_num4;
     
    	hpid = hp_nestid(b[2], b[3], 13);
       // htmid =lookupID(b[2], b[3], 13)/65536;
        // 16  64  256 1024  4096 16384 65536 
        // 11  10    9    8     7    6     5
    	q_num = hpid/256;
        
//        q_num = htmid/65536;
//        double ra_plus = b[2]+R_A;
//      	double ra_minus = b[2]-R_A;
//      	double dec_plus = b[3]+R_B;
//      	double dec_minus = b[3]-R_B;
//      	int hpid1 = hp_nestid(ra_plus, b[3], 13);
//      	int hpid2 = hp_nestid(ra_minus, b[3], 13);
//      	int hpid3 = hp_nestid(b[2], dec_plus, 13);
//      	int hpid4 = hp_nestid(b[2], dec_minus, 13);  
//        q_num1 = hpid1/65536;
//        q_num2 = hpid2/65536;
//        q_num3 = hpid3/65536;
//        q_num4 = hpid4/65536;
        
        
        //not have the key in map
	  if(recordmap2.find(q_num)==recordmap2.end()){ 
          	vector<vector<double>> temp;
   		temp.push_back(b);
         	recordmap2.insert(pair<int,vector<vector<double>>>(q_num, temp));
         	temp.clear();
     	  }	        
     	//already have the key in map
     	  else{
          	recordmap2[q_num].push_back(b);
     	  }
        //---------------------------------------------------
//        if (q_num1!=q_num){
//           if(recordmap2.find(q_num1)==recordmap2.end()){ 
//     		    vector<vector<double>> temp;
//     		    temp.push_back(b);
//            recordmap2.insert(pair<int,vector<vector<double>>>(q_num1, temp));
//            temp.clear();
//       		}	        
//       		//already have the key in map
//       		else{
//             recordmap2[q_num1].push_back(b);
//       		}
//        }
//        
//        if (q_num2 != q_num){
//           if(recordmap2.find(q_num2)==recordmap2.end()){ 
//     		    vector<vector<double>> temp;
//     		    temp.push_back(b);
//            recordmap2.insert(pair<int,vector<vector<double>>>(q_num2, temp));
//            temp.clear();
//       		}	        
//       		//already have the key in map
//       		else{
//             recordmap2[q_num2].push_back(b);
//       		}
//        }
//         
//         if (q_num3 != q_num){
//           if(recordmap2.find(q_num3)==recordmap2.end()){ 
//     		    vector<vector<double>> temp;
//     		    temp.push_back(b);
//            recordmap2.insert(pair<int,vector<vector<double>>>(q_num3, temp));
//            temp.clear();
//       		}	        
//       		//already have the key in map
//       		else{
//             recordmap2[q_num3].push_back(b);
//       		}
//         }
//          if (q_num4 != q_num){
//           if(recordmap2.find(q_num4)==recordmap2.end()){ 
//     		    vector<vector<double>> temp;
//     		    temp.push_back(b);
//            recordmap2.insert(pair<int,vector<vector<double>>>(q_num4, temp));
//            temp.clear();
//       		}	        
//       		//already have the key in map
//       		else{
//             recordmap2[q_num4].push_back(b);
//       		}
//        }
       //---------------------------------------------------
        b.clear(); 
   
    }  
    in.close();
	return recordmap2;
}

//get the index list of the divided blocks
vector<vector<int>> get_idlist(unordered_map<int,vector<vector<double>>> recordmap,vector<vector<int>> idlist){
	unordered_map<int,vector<vector<double>>>::iterator iter;
  	#pragma omp parallel for
  	for (iter=recordmap.begin();iter!=recordmap.end();iter++){
  		vector<int> temp;
  		temp.push_back(iter->first);
  		temp.push_back(iter->second.size());
		idlist.push_back(temp);
  		temp.clear();
  	}
  //sort(idlist.begin(), idlist.end());
	return idlist;
}

//get the shared index list of divided blocks of record A and B
vector<int> get_shared_id(unordered_map<int,vector<vector<double>>> recordmapA, unordered_map<int,vector<vector<double>>> recordmapB){
    vector<int> sharedlist;
    unordered_map<int,vector<vector<double>>>::iterator iterA;
    unordered_map<int,vector<vector<double>>>::iterator iterB;
    #pragma omp parallel for
    for (iterA=recordmapA.begin();iterA!=recordmapA.end();iterA++){
        if ((iterB = recordmapB.find(iterA->first)) != recordmapB.end()){
           sharedlist.push_back(iterA->first);
        }
    }
    //sort(sharedlist.begin(), sharedlist.end());
    return sharedlist;
}


vector<size_t> get_cudamalloc_size(vector<vector<unsigned int>> offset_band){
  set<unsigned int> tmp1;
  set<unsigned int> tmp2;
  set<unsigned int> tmp3;
  set<unsigned int> tmp4;
  
  vector<size_t> maxband;
  vector<size_t> maxnbytes;
  #pragma omp parallel for
  for (int m=0; m<offset_band.size(); m++){
      unsigned int data_x_band=offset_band[m][1];
      tmp1.insert(data_x_band);
      unsigned int data_y_band=offset_band[m][2];
      tmp2.insert(data_y_band);
     
  }
  maxband.push_back(*tmp1.rbegin());
  maxband.push_back(*tmp2.rbegin());

  size_t nBytes = (maxband[0])*sizeof(double);
  size_t nBytes2 = (maxband[1])*sizeof(double);
  size_t nBytes3 = (maxband[0])*(maxband[1])*sizeof(int);
  size_t nBytes4 = (maxband[0])*sizeof(int*);
  size_t nBytes5 = (maxband[1])*sizeof(int*);
  maxnbytes.push_back(nBytes);
  maxnbytes.push_back(nBytes2);
  maxnbytes.push_back(nBytes3);
  maxnbytes.push_back(nBytes4);
  maxnbytes.push_back(nBytes5);
  return maxnbytes;
}

vector<vector<double>> result_output(unsigned int data_x_band,unsigned int data_y_band,int  *h_out_dis,
                                    double *h_in_ra_x,double *h_in_dec_x, double*h_in_ra_y, double*h_in_dec_y,
                                    unsigned int data_x_offset,unsigned int data_y_offset,vector<vector<double>> matchresult){
    vector<double> tempresult;
    #pragma omp parallel for
    for (unsigned int ty = 0; ty<data_y_band; ty++){
  		for (unsigned int tx = 0; tx<data_x_band; tx++){
  			unsigned int tidx = (ty)*(data_x_band)+(tx);
                        if (h_out_dis[tidx]==1){
				double rax=h_in_ra_x[tx+data_x_offset];
				double decx=h_in_dec_x[tx+data_x_offset];
				double ray=h_in_ra_y[ty+data_y_offset];
				double decy=h_in_dec_y[ty+data_y_offset];
				tempresult.push_back(rax);
				tempresult.push_back(decx);
				tempresult.push_back(ray);
				tempresult.push_back(decy);
				tempresult.push_back(distance(rax,decx,ray,decy));
				matchresult.push_back(tempresult);
				tempresult.clear();
              	 	}
            	}
       }
     return matchresult;
}






