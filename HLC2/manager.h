#include "helper.h"
#include <omp.h>



int hp_nestid(double ra, double dec, int order);

unordered_map<int,vector<vector<double>>> read_to_unordered(string path,unordered_map<int,vector<vector<double>>> recordmap2);

vector<vector<int>> get_idlist(unordered_map<int,vector<vector<double>>> recordmap,vector<vector<int>> idlist);

vector<int> get_shared_id(unordered_map<int,vector<vector<double>>> recordmapA, unordered_map<int,vector<vector<double>>> recordmapB);

vector<size_t> get_cudamalloc_size(vector<vector<unsigned int>> offset_band);
                    
vector<double> parameter_decided();

int result_output(unsigned int data_x_band,unsigned int data_y_band,int *h_out_dis, double *h_in_ra_x,double *h_in_dec_x, double*h_in_ra_y, double*h_in_dec_y,unsigned int data_x_offset,unsigned int data_y_offset,FILE *fp3);