/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/* chealpix.h
 *
 */

#ifndef __CHEALPIX_H__
#define __CHEALPIX_H__

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* --------------------- */
/* Function Declarations */
/* --------------------- */

/* pixel operations */
/* ---------------- */
void ang2pix_nest(const long nside, double theta, double phi, long *ipix);
void ang2pix_ring(const long nside, double theta, double phi, long *ipix);

void mk_pix2xy(int *pix2x, int *pix2y);
void mk_xy2pix(int *x2pix, int *y2pix);


#endif /* __CHEALPIX_H__ */

