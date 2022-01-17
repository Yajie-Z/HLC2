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


double distance(double rax,double decx,double ray, double decy){
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



