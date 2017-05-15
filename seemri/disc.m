function iv = disc(R, res, T1, T2)
% DISC - Create disc phantom for SeeMRI
%   IV = DISC(R, RES) creates a imaging volume width a disc
%   of radius R and RES resolution.
%
%   IV = DISC(R, RES, T1, T2) also set T1 and T2 time constants.
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

if nargin < 3
    T1 = 0.8;
end
if nargin < 4
    T2 = 0.07;
end

u = @(x,y,R) double(sqrt(x.^2+y.^2)<=R);

x = -R:res:R;
y = -R:res:R;

[xs, ys] = meshgrid(x, y);

Mz0 = u(xs, ys, R);
Mz0(ys>0) = Mz0(ys>0)*0.5;

iv = ImagingVolume(x, y, T1, T2, Mz0);


