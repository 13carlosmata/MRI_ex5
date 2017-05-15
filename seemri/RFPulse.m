%RFPULSE - Create SeeMRI RF pulse object
%
%   RF = RFPulse(EF, A, F, PHI, TP) creates an RF pulse with envelope
%   function EF, amplitude A, frequency F (in Hz), phase PHI (in
%   radians) and pulse time TP (in seconds).
%
%   RF = RFPulse(EF, A, F, PHI, TP, T0) adds pulse start time T0 
%   (in seconds, default 0).
%
%   RF = RFPulse(EF, A, F, PHI, TP, T0, DT) adds simulation time step DT.
%   The time step should be small enough that the amplitude of the
%   envelope is approximately constant durting the time interval.
%
%   The RFPulse object represents an RF pulse according to: 
%
%      B1(t) = EF(t).*(t >= T0).*(t <= T0+TP) ...
%              .*(cos(2*pi*F*t+PHI)*x_vec - sin(2*pi*F*t+PHI)*y_vec)
%   
%    or in complex notation:
%
%      B1(t) = EF(t).*(t >= T0).*(t <= T0+TP).*exp(-i*(2*pi*F+PHI))
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

classdef RFPulse
    properties
        EnvelopeFunction
        Magnitude
        Frequency
        Phases
        Duration
        StartTime = 0;
        TimeStep
    end
    methods
        function rf = RFPulse(ef, A, f, phi, tp, t0, dt)
            if nargin > 0
                rf.EnvelopeFunction = ef;
                if nargin <= 5
                    t0 = 0;
                end             
                rf.Magnitude = A.*ones(size(t0));
                rf.Frequency = f;
                rf.Phases = phi.*ones(size(t0));
                rf.Duration = tp.*ones(size(t0));                
                rf.StartTime = t0;
                if nargin > 6
                    rf.TimeStep = dt;
                else
                    rf.TimeStep = rf.Duration/50;
                end
            end
        end
        function A = Amplitude(rf, t)
            A = 0;
            for ind = 1:length(rf.StartTime)
                A = A + rf.Magnitude(ind)...
                    .*rf.EnvelopeFunction(t-rf.StartTime(ind))...
                    .*(t>=rf.StartTime(ind)) ...
                    .*(t<=(rf.StartTime(ind)+rf.Duration(ind)));
            end
        end
        function phi = Phase(rf, t)        
            phi = rf.Phases(find(t >= rf.StartTime, 1, 'last'));
        end            
        function plot(rf, t1, t2)
        %PLOT - Plot envelope function
        %   PLOT(RF) plots the envelope function.
        %   
        %   PLOT(RF, T1, T2) plots the envelope function betweem times T1
        %   and T2.            
            if nargin < 2
                t1 = min(rf.StartTime);
                t2 = max(rf.StartTime + 2*rf.Duration);
            end
            ts = t1:(t2-t1)/200:t2;
            for ind = 1:length(rf.StartTime)
                ts = [ts ...
                      rf.StartTime(ind):rf.TimeStep:rf.StartTime(ind)+...
                      rf.Duration(ind)];
            end
            vals = [zeros(1, length(rf.StartTime)) rf.Amplitude(ts) ...
                   zeros(1, length(rf.StartTime))];
            ts = [rf.StartTime ts rf.StartTime+rf.Duration];
            [ts, sind] = sort(ts);
            vals = vals(sind);
            ind = find(ts >= t1 & ts <= t2);
            ts = ts(ind);
            vals = vals(ind);
            plot(ts, vals);
        end
        function h = mysubplot(rf, t1, t2, mn, mx)
            if nargin < 2
                t1 = min(rf.StartTime);
                t2 = max(rf.StartTime + 2*rf.Duration);
            end
            ts = t1:(t2-t1)/200:t2;
            for ind = 1:length(rf.StartTime)
                ts = [ts ...
                      rf.StartTime(ind):rf.TimeStep:rf.StartTime(ind)+...
                      rf.Duration(ind)];
            end
            vals = [zeros(1, length(rf.StartTime)) rf.Amplitude(ts) ...
                   zeros(1, length(rf.StartTime))];
            ts = [rf.StartTime ts rf.StartTime+rf.Duration];
            [ts, sind] = sort(ts);
            vals = vals(sind);
            ind = find(ts >= t1 & ts <= t2);
            ts = ts(ind);
            vals = vals(ind);
            if nargin > 3
                sigrange = (max(vals)-min(vals));
                if sigrange ~= 0
                    vals = (vals - min(vals))/sigrange*(mx-mn)+mn;
                else
                    vals = (vals - min(vals))+(mn+mx)/2;
                end
            end
            h = plot(ts, vals, 'k');                
        end
        function plot2(rf, t1, t2)
            for ind = 1:length(rf.StartTime)
                if nargin < 2
                    t1 = rf.StartTime(ind);
                    t2 = rf.StartTime(ind) + 2*rf.Duration(ind);
                end
                ts = t1:rf.TimeStep:t2;
                vals = rf.Amplitude(ts);
                cm = colormap('jet');
                patch([ts ts(end:-1:1)], [vals -vals(end:-1:1)], ...
                      cm(floor(rf.Phases(ind)/(2*pi/size(cm,1)))+1,:), ...
                      'EdgeColor', 'k');
            end
        end
        function [Y, fs] = powspec(rf)
        %POWSPEC - Plot power spectrum of envelope function 
        %   POWSPEC(RF) plots the power spectrum of the envelope function
        %   of RF.            
        %
        %   [Y, FS] = POWSPEC(RF) plots and returns the fourier transform
        %   Y at frequencies FS.            
            t1 = rf.StartTime - 10*rf.Duration;
            t2 = rf.StartTime + 10*rf.Duration;
            ts = t1:(t2-t1)/499:t2;
            fs = 1/(ts(2)-ts(1))*linspace(-.5,.5,length(ts));
            fs = fs-(fs(2)-fs(1))/2;
            vals = rf.EnvelopeFunction(ts).*(ts>=0)...
                   .*(ts<=(rf.Duration(1)));
            Y = fftshift(fft(vals));
            plot(fs, abs(Y));
            xlabel('Frequency (Hz)')
            ylabel('|F\{B_1^e\}|')
            %powspec(rf.EnvelopeFunction(ts-rf.StartTime), 1/(ts(2)-ts(1)));
        end
    end
end

