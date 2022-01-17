#include "helper.h"
#include <omp.h>



int hp_nestid(double ra, double dec, int order);

unordered_map<int,vector<vector<double>>> read_to_unordered(string path,unordered_map<int,vector<vector<double>>> recordmap2);

vector<vector<int>> combineblocks(int threshold,vector <vector<int>> combine,vector<int> sharedlist, vector<vector<int>> idlistA, vector<vector<int>>                                        idlistB,int p,int n,int m, vector <int> temp,int result);

vector<vector<int>> get_idlist(unordered_map<int,vector<vector<double>>> recordmap,vector<vector<int>> idlist);

vector<int> get_shared_id(unordered_map<int,vector<vector<double>>> recordmapA, unordered_map<int,vector<vector<double>>> recordmapB);

vector<size_t> get_cudamalloc_size(vector<vector<unsigned int>> offset_band);
                    
vector<double> parameter_decided();

     

