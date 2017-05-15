function iv = waterandfat(R, res)
% WATERANDFAT - Create water and fat phantom for imagin volume
%   IV = WATERANDFAT(R, RES) creates a disc, half water, half fat,
%   of radius R and resolution (pixel size in mm) RES.
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

u = @(x,y,R) double(sqrt(x.^2+y.^2)<=R);

x = -R:res:R;
y = -R:res:R;

[xs, ys] = meshgrid(x, y);

Mz0 = u(xs, ys, R);
Mz0(xs<0) = Mz0(xs<0)*0.5;

T1 = (xs<0)*2.3 + (xs>=0)*0.25;
T2 = (xs<0)*0.6 + (xs>=0)*0.07;

iv = ImagingVolume(x, y, T1, T2, Mz0, 'dB0', (xs>=0)*-3.35e-6);

