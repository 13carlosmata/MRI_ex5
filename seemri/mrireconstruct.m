function [imout, x, y] = mrireconstruct(S, kmax, varargin)
% MRIRECON - Reconstruct MRI image from k-space signal
%   [IM, X, Y] = MRIRECONSTRUCT(S, KMAX, ...) returns the
%   reconstructed image IM from the MRI signal S measured to a
%   maximum k-value of KMAX. X and Y contains the pixel
%   coorrdinates of the reconstructed image.
%
%   Options:
%     'Filter'    Set filter window to of: 'Hamming' (default), 
%                 'Gauss', 'none'.
%     'Plot'      Plot image: true or false (default)
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

ip = inputParser;
ip.FunctionName = 'mrireconstruct';
addOptional(ip, 'Filter', 'Hamming');
addOptional(ip, 'Plot', false);
ip.parse(varargin{:});

[M,N] = size(S);

switch ip.Results.Filter
  case 'Hamming'
    Wk = 0.54 + 0.46*cos(pi*(-N/2:N/2-1)/N);
  case  'Gauss'
    sigma = 0.4;
    Wk = exp(-((-N/2:N/2-1)/N).^2/sigma^2);
  case 'none'
    Wk = ones(1,N);
  otherwise
    disp('mrireconstruct: Unrecognized window, using Hamming window.')
    Wk = 0.54 + 0.46*cos(pi*(-N/2:N/2-1)/N);
end

x = (-N/2:1:N/2-1)/kmax/2;
y = (-M/2:1:M/2-1)/kmax/2;

im = abs(ifft2(fftshift(S.*(Wk'*Wk)...
                        .*((-1).^((1:M)'*ones(1,N) ...
                                  +ones(M,1)*(1:N))))));
if ip.Results.Plot
    imagesc(x, y, im)
    colormap gray
    axis image
end

if nargout > 0
    imout = im;
end