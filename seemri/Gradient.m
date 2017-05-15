%GRADIENT - Construct SeeMRI Gradient object
%
%    G = GRADIENT(T, A) creates a gradient sequence with amplitude
%    A(ind) at T(ind) <= t < T(ind+1) for ind = 1,... , end-1, and
%    A(end) for t >= T(end).
%
%    If A is a cell array then A{ind} can either be an array with the
%    gradient amplitude for each repetition or a function of the
%    repetion number N.
%  
%    Examples 
%       g1 = Gradient([0 0.1], [30e-3 0]);
%       g2 = Gradient([0 0.1], {[-30e-3:10e-3:20e-3], 0});
%       g3 = Gradient([0 0.1], {@(n) (n-4)*30e-3/3, 0});
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

classdef Gradient
    properties
        Times
        Amplitudes
    end
    methods
        function g = Gradient(ts, ampls)
            if nargin > 0
                g.Times = ts;
                g.Amplitudes = ampls;
            end
        end
        function G = Amplitude(g, t, n)
        % AMPLITUDE - Return amplitude at time 
        %   A = AMPLITUDE(G, T) returns the amplitude of G at time T.
        %   A = AMPLITUDE(G, T, N) returns the amplitude of G at
        %   time T for repetition N.            
            if nargin < 3
                n = 1;
            end
            if t < g.Times(1)
                G = 0;
            else 
                if iscell(g.Amplitudes)
                    a = g.Amplitudes{find(t >= g.Times, 1, 'last')};
                    if isnumeric(a)
                        if length(a) > 1
                            G = a(n);
                        else
                            G = a;
                        end
                    else
                        G = a(n);
                    end
                else
                    G = g.Amplitudes(find(t >= g.Times, 1, 'last'));
                end
            end
        end        
        function h = plot(g, t1, t2, n)
        %PLOT - Plot gradient
        %  PLOT(G) plots gradient object G.
        %  PLOT(G, T1, T2) plots G between times T1 and T2.
        %  PLOT(G, T1, T2, N) plots G between times T1 and T2 for 
        %  repetition N.
            if nargin < 2
                t1 = 0;
                t2 = 2*g.Times(end);
            end
            if nargin < 4
                n = 1;
            end
            ts = sort([t1 t2 g.Times-1e-9 g.Times+1e-9]);
            ts = ts(ts >= t1 & ts <= t2);
            h = plot(ts, arrayfun(@(t) Amplitude(g,t,n), ts));
        end
%         function subplot(g, t1, t2, mn, mx, n)
%             if nargin < 2
%                 t1 = 0;
%                 t2 = 2*g.Times(end);
%             end
%             if nargin < 6
%                 n = 1
%             end
%             ts = sort([t1 t2 g.Times-1e-9 g.Times+1e-9]);
%             ts = ts(ts >= t1 & ts <= t2);
%             vals = arrayfun(@(t) Amplitude(g,t,n), ts);
%             range = (max(vals)-min(vals));
%             if range ~= 0
%                 vals = (vals - min(vals))/range*(mx-mn)+mn;
%             else 
%                 vals = zeros(size(ts))+(mx+mn)/2;
%             end
%             plot(ts, vals, 'k');
%         end
        function h = subplot2(g, t1, t2, level, scale, n)
            if nargin < 2
                t1 = 0;
                t2 = 2*g.Times(end);
            end
            if nargin < 6
                n = 1;
            end
            ts = sort([t1 t2 g.Times-1e-9 g.Times+1e-9]);
            ts = ts(ts >= t1 & ts <= t2);
            vals = arrayfun(@(t) Amplitude(g,t,n), ts);
            h = plot(ts, scale*vals+level, 'k');
        end
        function [Gmin, Gmax] = AmplitudeRange(g, t1, t2, n)
            if nargin < 2
                t1 = 0;
                t2 = 2*g.Times(end);
            end
            if nargin < 4
                n = 1;
            end
            ts = sort([t1 t2 g.Times-1e-9 g.Times+1e-9]);
            ts = ts(ts >= t1 & ts <= t2);
            vals = arrayfun(@(t) Amplitude(g,t,n), ts);
            Gmin = min(vals);
            Gmax = max(vals);
        end
        function [Kmin, Kmax] = KRange(g, n)
            if nargin < 2
                n = 1;
            end
            ts = sort([g.Times-1e-9 g.Times+1e-9]);
            vals = arrayfun(@(t) Amplitude(g,t,n), ts);
            kvals = cumsum((vals(1:end-1)+vals(2:end))/2.*diff(ts));
            Kmin = min(kvals);
            Kmax = max(kvals);
        end        
    end
end