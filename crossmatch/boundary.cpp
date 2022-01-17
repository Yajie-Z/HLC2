#include "helper.h"
#include "constant.h"




//record time
double cpuSecond(){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

void file_to_string(vector<string> &record, const string& line, char delimiter)
{
    int linepos=0;
    char c;
    int linemax=line.length();
    string curstring;
    record.clear();

    while(linepos<linemax)
    {
        c = line[linepos];
        if(isdigit(c)||c=='.'){
            curstring+=c;
        }
        else if(c==delimiter&&curstring.size()){
            record.push_back(curstring);
            curstring="";
        }
        ++linepos;
    }
    if(curstring.size())
        record.push_back(curstring);
    return;
}

float string_to_float(string str){
    int i=0,len=str.length();
    float sum=0;
    while(i<len){
        if(str[i]=='.') break;
        sum=sum*10+str[i]-'0';
        ++i;
    }
    ++i;
    float t=1,d=1;
    while(i<len){
        d*=0.1;
        t=str[i]-'0';
        sum+=t*d;
        ++i;
    }
    return sum;
}


// read the star records from csv file to [format: unordered_map<hp index, corresponding star records> ]
vector<vector<double>> read_to_unordered(string path,vector<vector<double>> cataloglist){
	  vector<double>b;
    vector<string> row;
    
    string line;
    string filename;
    ifstream in(path);
    if (in.fail())  { cout << "File not found" <<endl; }
    
    while(getline(in, line)  && in.good() )
    {
        file_to_string(row, line, ',');  //把line里的单元格数字字符提取出来，“,”为单元格分隔符
        #pragma omp parallel for
        for(int i=0, leng=row.size(); i<leng; i++){
              b.push_back(string_to_float(row[i]));
        }
        cataloglist.push_back(b);
  
        b.clear(); 
   
    }  
    in.close();
	return cataloglist;
}

double pingfang(double a)
{
	return a*a;
}

vector<vector<double>> cataloglistA;
vector<vector<double>> cataloglistB;
int resultnum=0;

int main(int argc, char **argv){

    double AllStart = cpuSecond();
    double iStart, iElaps;
    double distance;
    double ra1,ra2,dec1,dec2;
    cataloglistA = read_to_unordered("data/sdsstest.csv",cataloglistA);
    iElaps = cpuSecond() - AllStart;
	  printf("A read to unorded time is %f s\n", iElaps);
    cataloglistB = read_to_unordered("data/twomasstest.csv",cataloglistB);
    iElaps = cpuSecond() - AllStart;
	  printf("A+B read to unorded time is %f s\n", iElaps);
    printf("A linenum: %d\n",cataloglistA.size());
    printf("B linenum: %d\n",cataloglistB.size());
    
    for(int i=0; i<cataloglistA.size(); i++){
      for(int j=0; j<cataloglistB.size(); j++){
         ra1 = cataloglistA[i][2];
         ra2 = cataloglistB[j][2];
         dec1 = cataloglistA[i][3];
         dec2 = cataloglistB[j][3];
        
        distance = pingfang(((double)(ra1-ra2))*((double)cos((dec1+dec2)/2.0/180.0*PI))) + pingfang((double)(dec1-dec2));
        if( distance < 9.0 * ( R_A*R_A + R_B*R_B ) ){
           resultnum++;
        }
      
      }
      printf("%d",resultnum);
    }
    
    
    
     
  	iElaps = cpuSecond() - AllStart;
    printf("result: %d\n", resultnum);
  	printf("[Info]All time is %f s\n", iElaps);	
   
  	return 0;
}