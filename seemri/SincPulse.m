%SINCPULSE - Create a SeeMRI sinc RF pulse object
%   
%   RF = SINCPULSE(B1, F, PHI, TP) creates a sinc RF pulse with
%   amplitude B1, frequency F, phase PHI and pulse time TP. The
%   envelope function is normalized such that it integrates to B1*TP.
%
%   RF = SINCPULSE(B1, F, PHI, TP, T0) adds pulse start time T0.
%
%   See also RFPULSE, RFPULSE/PLOT, RFPULSE/POWSPEC
%
%   Copyright (c) 2010 Peter Nillius (nillius@mi.physics.kth.se)
%
 
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the
%   License. You may obtain a copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
%   implied. See the License for the specific language governing
%   permissions and limitations under the License.

classdef SincPulse < RFPulse
    properties (SetAccess = private)
        BandWidth
    end
    methods
       function rf = SincPulse(B1, f, phi, tp, t0)
           if nargin == 0
               EnvelopeAmplitude = 0;
           else
               EnvelopeAmplitude = B1;
           end
           if nargin < 5
               t0 = 0;
           end
           nzeros = 6;
           df = 2*nzeros/tp;
           A = tp/quad(@(t) 2*sinc(df*t), 0, tp/2, 1e-12);
           rf = rf@RFPulse(@(t) A*sinc(df*(t-tp/2)), B1, f, phi, tp, t0, ...
                           tp/200);
           rf.BandWidth = df;
       end
   end
end
