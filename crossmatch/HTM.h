/*
 * HTM.h
 *
 *  Created on: 2009-11-7
 *     
 */

#ifndef HTM_H_
#define HTM_H_

#define Pi 		3.1415926535897932384626433832795
#define Pr 		3.1415926535897932384626433832795/180.0
#define Epsilon 	1.0E-15
#define gEpsilon 	-1.0E-15
#define sqrt3		1.7320508075688772935
#define IDSIZE		64
#define	HTMNAMEMAX	32

#define iS2		0
#define iN1		1
#define iS1		2
#define iN2		3
#define iS3		4
#define iN0		5
#define iS0		6
#define iN3		7

struct Base
{
	char name[3];
	int ID;
	int v1, v2, v3;
};

int startpane(double* v1, double* v2, double* v3,
	double xin, double yin, double zin,
	char* name);

void m4_midpoint(double* v1, double* v2, double* w);
char* append(char* str, char c);
void handleError(const char* error);
void copy_vec(double* d, double* s);
int isinside(double* p, double* v1, double* v2, double* v3);
double* radecToVector(double ra, double dec, double* vec);

int lookup(double x, double y, double z, int depth, char* name);
long nameToId(char* name);

int lookupID(double ra, double dec, int depth);
char* lookupName(double ra, double dec, int depth, char* name);


#endif /* HTM_H_ */
