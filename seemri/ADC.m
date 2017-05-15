% ADC - Create a SeeMRI analog-to-digital converter object
%
%   A = ADC(TEND, DT) creates an object which samples from time 0 to
%   TEND with sampling interval DT.
%  
%   A = ADC(TSTART, TEND, DT) creates an object which samples from
%   time TSTART to TEND with sampling interval DT. The exact sampling
%   times are adapted such that the number of sampling points, N, are
%   even and point number N/2+1 is in the center of the interval. This
%   is suitable for reconstruction when sampling around an echo.
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

classdef ADC
    properties
        Times
    end
    methods
        function adc = ADC(varargin)
            switch length(varargin)
              case 2
                t_end = varargin{1};
                dt = varargin{2};
                adc.Times = 0:dt:t_end;
              case 3                
                t_start = varargin{1};
                t_end = varargin{2};
                dt = varargin{3};
                N = floor((t_end-t_start)/dt) + 1;
                % make N even
                N = N - mod(N,2); 
                % Let sample number N/2+1 be at center (echo time)
                t_start = (t_start+t_end)/2 - N/2*dt;
                t_end = t_start + (N-1)*dt;
                adc.Times = t_start:dt:t_end;
            end
        end
        function h = subplot(adc, t1, t2, mn, mx)
            if nargin < 2
                t1 = 0;
                t2 = 2*adc.Times(end);
            end
            ts = [t1 adc.Times(1) adc.Times adc.Times(end) t2];
            vals = [0 0 ones(1, length(adc.Times)) 0 0];
            range = (max(vals)-min(vals));
            if range ~= 0
                vals = (vals - min(vals))/range*(mx-mn)+mn;
            else 
                vals = zeros(size(ts))+(mx+mn)/2;
            end
            h = plot(ts, vals, 'k');
        end

    end
end