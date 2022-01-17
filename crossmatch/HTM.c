/*
 * HTM.c
 *
 *  Created on: 2009-11-7
 *
 */

#include "HTM.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
//#include <stdlib.h>

const int N_indexes[4][3] = {
	{1, 0, 4},	// N0
	{4, 0, 3},	// N1
	{3, 0, 2},	// N2
	{2, 0, 1}	// N3
};

const int S_indexes[4][3] = {
	{1, 5, 2},	// S0
	{2, 5, 3},	// S1
	{3, 5, 4},	// S2
	{4, 5, 1},	// S3
};

const int anchor[6][3] = {
	{0,  0,  1},	// 0	a0
	{1,  0,  0},	// 1	a1
	{0,  1,  0},	// 2	a2
	{-1, 0,  0},	// 3	a3
	{0, -1,  0},	// 4	a4
	{0,  0, -1}	// 5	a5
};

const struct Base bases[] = {
	{"S2", 10, 3, 5, 4},
	{"N1", 13, 4, 0, 3},
	{"S1", 9, 2, 5, 3},
	{"N2", 14, 3, 0, 2},
	{"S3", 11, 4,5,1},
	{"N0", 12, 1, 0, 4},
	{"S0", 8, 1, 5, 2},
	{"N3", 15, 2, 0, 1}
};


int startpane(double* v1, double* v2, double* v3,
	double xin, double yin, double zin,
	char* name)
{
	//double tvec[3];
	int baseID;
	int baseindex = 0;
	if ((xin > 0) && (yin >= 0)){
	baseindex = (zin >= 0) ? iN3 : iS0;

	} else if ((xin <= 0) && (yin > 0)){
	baseindex = (zin >= 0) ? iN2 : iS1;

	} else if ((xin < 0) && (yin <= 0)){
	baseindex = (zin >= 0) ? iN1 : iS2;

	} else if ((xin >= 0) && (yin < 0)){
	baseindex = (zin >= 0) ? iN0 : iS3;
	} else {
	//ErrorHandler.handleError(ErrorHandler.PANIC);
	}


	baseID = bases[baseindex].ID;

	//tvec = anchor[bases[baseindex].v1];
	v1[0] = anchor[bases[baseindex].v1][0];
	v1[1] = anchor[bases[baseindex].v1][1];
	v1[2] = anchor[bases[baseindex].v1][2];

	//tvec = (double[])anchor[bases[baseindex].v2];
	v2[0] = anchor[bases[baseindex].v2][0];
	v2[1] = anchor[bases[baseindex].v2][1];
	v2[2] = anchor[bases[baseindex].v2][2];

	//tvec = (double[])anchor[bases[baseindex].v3];
	v3[0] = anchor[bases[baseindex].v3][0];
	v3[1] = anchor[bases[baseindex].v3][1];
	v3[2] = anchor[bases[baseindex].v3][2];

	//name.append(bases[baseindex].name);
	strcpy(name, bases[baseindex].name);

	return baseID;
}


void m4_midpoint(double* v1, double* v2, double* w)
{
	w[0] = v1[0] + v2[0];
	w[1] = v1[1] + v2[1];
	w[2] = v1[2] + v2[2];
	double tmp = sqrt(w[0] * w[0] + w[1] * w[1] + w[2]*w[2]);
	w[0] /= tmp;
	w[1] /= tmp;
	w[2] /= tmp;
}

char* append(char* str, char c)
{
	int n = strlen(str);
	str[n] = c;
	str[n+1] = 0;
	return str;
}

void copy_vec(double* d, double* s)
{
	int i;
	for(i=0; i<3; i++ )
	{
		d[i] = s[i];
	}
}

/**
 * for a given vector p is it contained in the triangle whose corners are
 *  given by the vectors v1, v2,v3.
 */
int isinside(double* p, double* v1, double* v2, double* v3) // p need not be nromalilzed!!!
{
	double crossp[3];


	crossp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	crossp[1] = v1[2] * v2[0] - v2[2] * v1[0];
	crossp[2] = v1[0] * v2[1] - v2[0] * v1[1];
	if (p[0] * crossp[0] + p[1] * crossp[1] + p[2] * crossp[2] < gEpsilon)
		return 0;

	crossp[0] = v2[1] * v3[2] - v3[1] * v2[2];
	crossp[1] = v2[2] * v3[0] - v3[2] * v2[0];
	crossp[2] = v2[0] * v3[1] - v3[0] * v2[1];
	if (p[0] * crossp[0] + p[1] * crossp[1] + p[2] * crossp[2] < gEpsilon)
		return 0;


	crossp[0] = v3[1] * v1[2] - v1[1] * v3[2];
	crossp[1] = v3[2] * v1[0] - v1[2] * v3[0];
	crossp[2] = v3[0] * v1[1] - v1[0] * v3[1];
	if (p[0] * crossp[0] + p[1] * crossp[1] + p[2] * crossp[2] < gEpsilon)
		return 0;

	return 1;
}


void handleError(const char* error)
{
	printf("error: %s\n", error);
}


double* radecToVector(double ra, double dec, double* vec)
{
	double v[3];
	double cd = cos( (dec * Pi) / 180.0);

	double diff;
	diff = 90.0 - dec;

	if (diff < Epsilon && diff > gEpsilon){
		v[0] = 1.0;
		v[1] = 0.0;
		v[2] = 1.0;
		copy_vec(vec, v);
		return vec;
	}

	diff = -90.0 - dec;
	if (diff < Epsilon && diff > gEpsilon){
		v[0] = 1.0;
		v[1] = 0.0;
		v[2] = -1.0;
		copy_vec(vec, v);
		return vec;
	}
	v[2] = sin((dec* Pi) / 180.0);
	double quadrant;
	double qint;
	int iint;
	quadrant = ra / 90.0; // how close is it to an integer?

	// if quadrant is (almost) an integer, force x, y to particular
	// values of quad:
	// quad,   (x,y)
	// 0       (1,0)
	// 1,      (0,1)
	// 2,      (-1,0)
	// 3,      (0,-1)
	// q>3, make q = q mod 4, and reduce to above
	// q mod 4 should be 0.

	qint = (int)(quadrant+0.5);
    double r = (double) (qint - quadrant);
    if (r<0){
        r=r*(-1);

    }
	if(r < Epsilon)
	{
		iint = (int) qint;
		iint %= 4;
		if (iint < 0) iint += 4;

		switch(iint)
		{
		  case 0:
			  v[0] = 1.0;
			  v[1] = 0.0;
			  break;
		  case 1:
			  v[0] = 0.0;
			  v[1] = 1.0;
			  break;
		  case 2:
			  v[0] = -1.0;
			  v[1] = 0.0;
			  break;
		  case 3:
			  v[0] = 0.0;
			  v[1] = -1.0;
			  break;
		  }
		copy_vec(vec, v);
		return vec;
	}

	v[0] = cos((ra * Pi) / 180.0) * cd;
	v[1] = sin((ra * Pi) / 180.0) * cd;
	copy_vec(vec, v);
	return vec;
}



int lookup(double x, double y, double z, int depth, char* name)
{

	int startID;

	char sname[80];

	double v1[3];
	double v2[3];
	double v0[3];
	double w1[3];
	double w2[3];
	double w0[3];
	double p[3];

	p[0] = x;
	p[1] = y;
	p[2] = z;

	// Get the ID of the level0 triangle, and its starting vertices
	startID = startpane(v0, v1, v2, x, y, z, sname);

	// Start searching for the children
	while(depth-- > 0)
	{
		m4_midpoint(v0, v1, w2);
		m4_midpoint(v1, v2, w0);
		m4_midpoint(v2, v0, w1);

		if (isinside(p, v0, w2, w1))
		{
			append(sname, '0');
			copy_vec(v1, w2);
			copy_vec(v2, w1);
			//printf("%s\n", sname);
		}
		else if (isinside(p, v1, w0, w2))
		{
			append(sname, '1');
			copy_vec(v0, v1);
			copy_vec(v1, w0);
			copy_vec(v2, w2);
		}
		else if (isinside(p, v2, w1, w0))
		{
			append(sname, '2');
			copy_vec(v0, v2);
			copy_vec(v1, w1);
			copy_vec(v2, w0);
		}
		else if (isinside(p, w0, w1, w2))
		{
			append(sname, '3');
			copy_vec(v0, w0);
			copy_vec(v1, w1);
			copy_vec(v2, w2);
		}
	}

	strcpy(name, sname);
	//return name.toString();
	return 0;
}


long nameToId(char* name)
{
	long out=0;
	int i;
	int siz = 0;

	if(name == 0 || strlen(name) == 0)			  // null pointer-name
		handleError("HTM name is null pointer-name");
	if(name[0] != 'N' && name[0] != 'S')  // invalid name
		handleError("HTM name is INVALIDNAME");

	siz = strlen(name);	   // determine string length
	// at least size-2 required, don't exceed max
	if(siz < 2)
		handleError("HTM name is INVALIDNAME");
	if(siz > HTMNAMEMAX)
		handleError("HTM name is INVALIDNAME");

	for(i = siz-1; i > 0; i--)
	{// set bits starting from the end
		if(name[i] > '3' || name[i] < '0') {// invalid name
			handleError("HTM name is INVALIDNAME");
	}
		out += ( (long)(name[i]-'0')) << 2*(siz - i -1);
	}

	i = 2;					 // set first pair of bits, first bit always set
	if(name[0]=='N') i++;	  // for north set second bit too
	long last = ((long)i << (2*siz - 2) );
	out += last;
	return out;
}

int lookupID(double ra, double dec, int depth)
{
	double vec[3];
	radecToVector(ra, dec, vec);
	//printf("%g,%g,%g\n", vec[0], vec[1], vec[2]);
	char name[80];
	lookup(vec[0], vec[1], vec[2], depth, name);

	return nameToId(name);
}

