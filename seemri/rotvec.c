/*
 * File:	$RCSfile$
 * Version:	$Revision$
 * Revised:	$Date$
 *
 * Created:	2010/04/30
 * Author:	Peter Nillius
 *
 * Description: 
 *
 * Copyright (c) 2010 Peter Nillius Licensed under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License
 * at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 *
 * ----------------------------------------------------------------------
 */

static char version[] = "@(#)$Id$\n\
@(#)$Copyright: (C) Copyright 2010 Peter Nillius $";

#include <mex.h>
#include <math.h>

void mexFunction(int noargs,              /* Antal utparametrar        */
		 mxArray *oargs[],        /* Ut-parameter-pekar-vektor */
		 int niargs,              /* Antal inparametrar        */
		 const mxArray *iargs[])  /* In-parameter-pekar-vektor */
{
    const int nargs = 3;
    int argind = 0;
    const mxArray *v, *z, *phi;
    double *vp, *zp, *phip, *vrotp;
    int c;
    double cosphi, sinphi, zdotv;

    if (niargs < nargs) {
	mexPrintf("need %d args.\n", nargs);
	return;	
    }
    
    v = iargs[argind++];
    z = iargs[argind++];
    phi = iargs[argind++];

    oargs[0] = mxCreateDoubleMatrix(mxGetM(v), mxGetN(v), mxREAL);

    vp = (double *)mxGetData(v);
    zp = (double *)mxGetData(z);
    phip = (double *)mxGetData(phi);
    vrotp = (double *)mxGetData(oargs[0]);
    for (c = 0; c < mxGetN(v); c++) {
        cosphi = cos(phip[c]);        
        sinphi = sin(phip[c]);
        zdotv = zp[0]*vp[0] + zp[1]*vp[1] + zp[2]*vp[2];
        vrotp[0] = vp[0]*cosphi + (zp[1]*vp[2]-zp[2]*vp[1])*sinphi 
            + zp[0]*zdotv*(1-cosphi);
        vrotp[1] = vp[1]*cosphi + (zp[2]*vp[0]-zp[0]*vp[2])*sinphi 
            + zp[1]*zdotv*(1-cosphi);
        vrotp[2] = vp[2]*cosphi + (zp[0]*vp[1]-zp[1]*vp[0])*sinphi 
            + zp[2]*zdotv*(1-cosphi);
        vp += 3;
        zp += 3;
        vrotp += 3;
    }

}



