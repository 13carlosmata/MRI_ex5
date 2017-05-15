%IMAGINGVOLUME - Create a SeeMRI imaging volume object
%
%   IV = IMAGINGVOLUME(X, Y, T1, T2) create and imaging volume in a
%   grid with coordinates X and Y on the x and y axis
%   respectively. The total number of magnetization vectors is equal
%   to LENGTH(X)*LENGTH(Y). T1 and T2 sets the relaxation constant at
%   each point and may be a single value or a matrix of size LENGTH(Y)
%   by LENGTH(X).
%
%   IV = IMAGINGVOLUME(X, Y, T1, T2, MZ0, ...) also sets the bulk
%   magnetization equilibrium value MZ0 (default 1).
%
%   Options:
%    'dB0Sigma'     Set the standard deviation of a normally
%                   distributed B0 inhomogeniety with normal The total
%                   field strength is calculated as B0*(1+dB0).
%
%    'dB0Gamma'     Set the gamma parmeter of the B0 inhomogeniety
%                   with a Cauchy-Lorentz distribution. The total field
%                   strength is calculated as B0*(1+dB0).
%
%    'dB0'          Set relative B0 deviation map. The total field
%                   strength is calculated as B0*(1+dB0).
%
%   See also IMAGINGVOLUME/TOEQUILIBRIUM, IMAGINGVOLUME/PLOT
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

classdef ImagingVolume < handle
    properties
        x
        y
        M
        Mz0
        T1
        T2
    end
    properties (SetAccess = private)
        dB0ratio
        mh
        ip
        plot_margin
        plot_M_scale
        fancyplot
    end
    methods
        function iv = ImagingVolume(x, y, T1, T2, Mz0, varargin)
            iv.ip = inputParser;
            iv.ip.FunctionName = 'ImagingVolume';
            addOptional(iv.ip, 'Annotation', []);
            addOptional(iv.ip, 'dB0Sigma', 0);        
            addOptional(iv.ip, 'dB0Gamma', 0);        
            addOptional(iv.ip, 'dB0', 0);        
            addOptional(iv.ip, 'M', []);     
            addOptional(iv.ip, 'PlotScale', 1);        
            iv.ip.parse(varargin{:});
            
            [xs, ys] = meshgrid(x, y);
            iv.x = xs(:)';
            iv.y = ys(:)';
            if nargin < 5
                iv.Mz0 = ones(size(iv.x));
            else
                iv.Mz0 = ones(size(iv.x)).*(Mz0(:)');
            end
            if length(iv.ip.Results.M) > 0
                iv.M = iv.ip.Results.M;
            else
                toEquilibrium(iv);
            end
            iv.T1 = ones(size(iv.x)).*T1(:)';
            iv.T2 = ones(size(iv.x)).*T2(:)';
            if iv.ip.Results.dB0Sigma == 0  
                % Cauchy-Lorentz distribution
                iv.dB0ratio = iv.ip.Results.dB0Gamma...
                    *tan(pi*unifrnd(0, 1, size(iv.x)));
            else
                % Normal distribution
                iv.dB0ratio = normrnd(0, iv.ip.Results.dB0Sigma, ...
                                      size(iv.x));
            end
            iv.dB0ratio = iv.dB0ratio + iv.ip.Results.dB0(:)';
            if length(x) == 1
                iv.plot_M_scale = 1/max(iv.Mz0);
                iv.plot_margin = 1;
            else
                iv.plot_margin = max((max(x)-min(x))/(length(x)-1), ...
                                     (max(y)-min(y))/(length(y)-1))/2;
                iv.plot_M_scale = iv.plot_margin/max(iv.Mz0);
            end
            iv.plot_M_scale = iv.plot_M_scale * iv.ip.Results.PlotScale;
        end
        function toEquilibrium(iv)
        %TOEQUILIBRIUM - Set magnetization to equilibrium
        %   TOEQUILIBRIUM(IV) sets the magnetization vectors in IV
        %   to their equilibrium values.
            iv.M = zeros(3,length(iv.x));            
            iv.M(3,:) = iv.Mz0;
        end
        function hout = plot(iv)
        %PLOT - Plot magnetization vectors of imaging volume
        %   PLOT(IV) plots the M vectors of IV.
            xmin = min(min(iv.x));
            xmax = min(max(iv.x));
            ymin = min(min(iv.y));
            ymax = min(max(iv.y));
            marg = iv.plot_margin;
            bh = plot3([[-marg -marg -marg -marg]+xmin ...
                        [marg marg marg marg]+xmax], ...
                       [-marg+ymin -marg+ymin marg+ymax marg+ymax ...
                        -marg+ymin -marg+ymin marg+ymax marg+ymax], ...
                       [-marg marg -marg marg -marg marg -marg marg], '.');
            hold on
            iv.fancyplot = false; %length(iv.x) == 1;
            if iv.fancyplot
                iv.mh = fancyvector([iv.x; iv.x+iv.plot_M_scale*iv.M(1,:)], ...
                                    [iv.y; iv.y+iv.plot_M_scale*iv.M(2,:)], ...
                                    [zeros(1, size(iv.M,2)); ...
                                    iv.plot_M_scale*iv.M(3,:)], ...
                                    0, max(iv.Mz0)/20);
                light
                lighting gouraud
                shading interp
                colormap([0 0 1])
            else
                iv.mh = plot3([iv.x; iv.x + iv.plot_M_scale*iv.M(1,:)], ...
                              [iv.y; iv.y + iv.plot_M_scale*iv.M(2,:)], ...
                              [zeros(1, size(iv.M,2)); ...
                               iv.plot_M_scale*iv.M(3,:)], 'b', ...
                              'LineWidth', 2);
            end
            h = iv.mh;
            
            if length(iv.ip.Results.Annotation) > 0
                if iscellstr(iv.ip.Results.Annotation)
                    text(iv.x, iv.y, zeros(1, size(iv.M,2))-0.1, ...
                         iv.ip.Results.Annotation, ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'top');
                else
                    text(iv.x, iv.y, zeros(1, size(iv.M,2))-0.1, ...
                         arrayfun(iv.ip.Results.Annotation, ...
                                  1:size(iv.M,2), iv.T1, iv.T2, iv.Mz0,...
                                  'UniformOutput', false), ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'top');
                end
            end
            hold off
            axis equal
            if nargout > 0
                hout = h;
            end            
        end
        function update_plot(iv, phi)
            if nargin < 2
                if iv.fancyplot
                    fancyvector([iv.x; ...
                                 iv.x+iv.plot_M_scale*iv.M(1,:)], ...
                                [iv.y; ...
                                 iv.y+iv.plot_M_scale*iv.M(2,:)], ...
                                [zeros(1, size(iv.M,2)); ...
                                 iv.plot_M_scale*iv.M(3,:)], ...
                                0, iv.Mz0(1)/20, iv.mh);
                else
                    arrayfun(@(h, Mx, My, Mz, x, y) ...
                             set(h, 'XData', [0 Mx]+x, ...
                                    'YData', [0 My]+y, ...
                                    'ZData', [0 Mz]), ...
                             iv.mh', iv.plot_M_scale*iv.M(1,:), ...
                             iv.plot_M_scale*iv.M(2,:), ...
                             iv.plot_M_scale*iv.M(3,:), iv.x, ...
                             iv.y);
                end
            else
                if iv.fancyplot
                    fancyvector([iv.x; ...
                                 (iv.x ...
                                  +cos(phi)*iv.plot_M_scale*iv.M(1,:) ...
                                  +sin(phi)*iv.plot_M_scale*iv.M(2,:))], ...
                                [iv.y; ...
                                 (iv.y ...
                                  -sin(phi)*iv.plot_M_scale*iv.M(1,:) ...
                                  +cos(phi)*iv.plot_M_scale*iv.M(2,:))], ...
                                [zeros(1, size(iv.M,2)); ...
                                 iv.plot_M_scale*iv.M(3,:)], ...
                                0, iv.Mz0(1)/20, iv.mh);
                else
                    arrayfun(@(h, Mx, My, Mz, x, y) ...
                             set(h, 'XData', [0 Mx*cos(phi)+My*sin(phi)]+x, ...
                                    'YData', [0 -Mx*sin(phi)+My*cos(phi)]+y, ...
                                    'ZData', [0 Mz]), ...
                             iv.mh', iv.plot_M_scale*iv.M(1,:), ...
                             iv.plot_M_scale*iv.M(2,:), ...
                             iv.plot_M_scale*iv.M(3,:), iv.x, iv.y);
                end
            end
        end        
    end
end

