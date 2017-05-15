/*
 * File:	$RCSfile$
 * Version:	$Revision$
 * Revised:	$Date$
 *
 * Created:	2010/07/15
 * Author:	Peter Nillius
 *
 * Description: 
 *
 * Copyright (c) 2010-2013 Peter Nillius 
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you
 * may not use this file except in compliance with the License. You
 * may obtain a copy of the License at
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
    const int nargs = 7;
    int argind = 0;
    const mxArray *v, *Beff, *T1, *T2, *Mz0;
    double *vp, *Beffp, *T1p, *T2p, *Mz0p, *nvp;
    float *vfp, *Befffp, *T1fp, *T2fp, *Mz0fp, *nvfp;
    int c;
    double cosphi, sinphi, zdotv, Beffmag, gamma, dt, dalpha, z[3], 
        R2, tmp1;
    float cosphif, sinphif, zdotvf, Beffmagf, gammaf, dtf, dalphaf, zf[3], 
        R2f, tmp1f;
    
    if (niargs < nargs) {
	mexPrintf("need %d args.\n", nargs);
	return;	
    }
    
    v = iargs[argind++];
    Beff = iargs[argind++];
    T1 = iargs[argind++];
    T2 = iargs[argind++];
    Mz0 = iargs[argind++];

    if (mxIsDouble(v) && mxIsDouble(Beff) && mxIsDouble(T1)
        && mxIsDouble(T2) && mxIsDouble(Mz0)) {    

        gamma = mxGetScalar(iargs[argind++]);
        dt = mxGetScalar(iargs[argind++]);
        
        oargs[0] = mxCreateDoubleMatrix(mxGetM(v), mxGetN(v), mxREAL);

        vp = (double *)mxGetData(v);
        Beffp = (double *)mxGetData(Beff);
        T1p = (double *)mxGetData(T1);
        T2p = (double *)mxGetData(T2);
        Mz0p = (double *)mxGetData(Mz0);
        nvp = (double *)mxGetData(oargs[0]);
        if (mxGetM(Beff) == 3) {
            for (c = 0; c < mxGetN(v); c++) {
                Beffmag = sqrt(Beffp[0]*Beffp[0] + Beffp[1]*Beffp[1] 
                               + Beffp[2]*Beffp[2]);
                if (Beffmag != 0) {
                    z[0] = Beffp[0]/Beffmag;
                    z[1] = Beffp[1]/Beffmag;
                    z[2] = Beffp[2]/Beffmag;
                    zdotv = z[0]*vp[0] + z[1]*vp[1] + z[2]*vp[2];
                
                    dalpha = -Beffmag*gamma*dt;
                    cosphi = cos(dalpha);        
                    sinphi = sin(dalpha);
                
                    tmp1 = zdotv*(1-cosphi);
                    nvp[0] = vp[0]*cosphi + (z[1]*vp[2]-z[2]*vp[1])*sinphi 
                        + z[0]*tmp1;
                    nvp[1] = vp[1]*cosphi + (z[2]*vp[0]-z[0]*vp[2])*sinphi 
                        + z[1]*tmp1;
                    nvp[2] = vp[2]*cosphi + (z[0]*vp[1]-z[1]*vp[0])*sinphi 
                        + z[2]*tmp1;
                } else {
                    nvp[0] = vp[0];
                    nvp[1] = vp[1];
                    nvp[2] = vp[2];
                }
            
                R2 = exp(-dt/T2p[c]);
                nvp[0] = nvp[0]*R2;
                nvp[1] = nvp[1]*R2;
                nvp[2] = Mz0p[c] - (Mz0p[c] - nvp[2])*exp(-dt/T1p[c]);
            
                vp += 3;
                Beffp += 3;
                nvp += 3;
            }
        } else {
            for (c = 0; c < mxGetN(v); c++) {
            
                dalpha = -Beffp[c]*gamma*dt;
                cosphi = cos(dalpha);        
                sinphi = sin(dalpha);
            
                R2 = exp(-dt/T2p[c]);
                nvp[0] = (cosphi*vp[0] - sinphi*vp[1])*R2;
                nvp[1] = (sinphi*vp[0] + cosphi*vp[1])*R2;
                nvp[2] = Mz0p[c] - (Mz0p[c] - vp[2])*exp(-dt/T1p[c]);
            
                vp += 3;
                nvp += 3;
            }
        }
    } else if (mxIsSingle(v) && mxIsSingle(Beff) && mxIsSingle(T1)
               && mxIsSingle(T2) && mxIsSingle(Mz0)) {    
        //mexPrintf("single precision.\n");

        gammaf = mxGetScalar(iargs[argind++]);
        dtf = mxGetScalar(iargs[argind++]);
        //mexPrintf("single precision. %g\n", gammaf);
        
        mwSize dims[] = {mxGetM(v), mxGetN(v)};
        oargs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
        
        vfp = (float *)mxGetData(v);
        Befffp = (float *)mxGetData(Beff);
        T1fp = (float *)mxGetData(T1);
        T2fp = (float *)mxGetData(T2);
        Mz0fp = (float *)mxGetData(Mz0);
        nvfp = (float *)mxGetData(oargs[0]);
        if (mxGetM(Beff) == 3) {
            for (c = 0; c < mxGetN(v); c++) {
                Beffmagf = sqrtf(Befffp[0]*Befffp[0] + Befffp[1]*Befffp[1] 
                                 + Befffp[2]*Befffp[2]);
                if (Beffmagf != 0) {
                    zf[0] = Befffp[0]/Beffmagf;
                    zf[1] = Befffp[1]/Beffmagf;
                    zf[2] = Befffp[2]/Beffmagf;
                    zdotv = zf[0]*vfp[0] + zf[1]*vfp[1] + zf[2]*vfp[2];
                
                    dalphaf = -Beffmagf*gammaf*dtf;
                    cosphif = cosf(dalphaf);        
                    sinphif = sinf(dalphaf);
                
                    tmp1f = zdotvf*(1-cosphif);
                    nvfp[0] = vfp[0]*cosphif 
                        + (zf[1]*vfp[2]-zf[2]*vfp[1])*sinphif + zf[0]*tmp1f;
                    nvfp[1] = vfp[1]*cosphif 
                        + (zf[2]*vfp[0]-zf[0]*vfp[2])*sinphif + zf[1]*tmp1f;
                    nvfp[2] = vfp[2]*cosphif 
                        + (zf[0]*vfp[1]-zf[1]*vfp[0])*sinphif + zf[2]*tmp1f;
                } else {
                    nvfp[0] = vfp[0];
                    nvfp[1] = vfp[1];
                    nvfp[2] = vfp[2];
                }
            
                R2 = expf(-dt/T2fp[c]);
                nvfp[0] = nvfp[0]*R2;
                nvfp[1] = nvfp[1]*R2;
                nvfp[2] = Mz0fp[c] - (Mz0fp[c] - nvfp[2])*expf(-dt/T1fp[c]);
            
                vfp += 3;
                Befffp += 3;
                nvfp += 3;
            }
        } else {

            for (c = 0; c < mxGetN(v); c++) {            
                dalphaf = -Befffp[c]*gammaf*dtf;
                cosphif = cosf(dalphaf);
                sinphif = sinf(dalphaf);
            
                R2f = expf(-dtf/T2fp[c]);
                nvfp[0] = (cosphif*vfp[0] - sinphif*vfp[1])*R2f;
                nvfp[1] = (sinphif*vfp[0] + cosphif*vfp[1])*R2f;
                nvfp[2] = Mz0fp[c] - (Mz0fp[c] - vfp[2])*expf(-dtf/T1fp[c]);
            
                vfp += 3;
                nvfp += 3;
            }
        }
        
    } else {
        mexPrintf("Error: mixed double and single precision.\n");
    }

}



