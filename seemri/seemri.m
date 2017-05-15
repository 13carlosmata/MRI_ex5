function [SN, tNs] = seemri(iv, B0, rf, varargin)
%SEEMRI - Apply pulse sequence and simulate Bloch equations.
%
%   [S, T] = SEEMRI(IV, B0, RF, GX, GY, ADC, ...) applies to the
%   imaging volume IV, at magnetic field strength B0, the pulse
%   sequence with pulse RF, x-gradient GX, y-gradient GY. An empty
%   array, [], can be supplied to either RF, GX or GY, if no pulse or
%   gradient is to be used.
%
%   [S, T] = SEEMRI(IV, B0, RF, GX, GY, ADC, TR, N, ...) applies the
%   pulse sequence N times at starting times 0, TR, 2*TR, ...,
%   (N-1)*TR.
%
%   [S, T] = SEEMRI(IV, B0, RF, ...) applies the pulse without
%   gradients and samples the signal from the beginning of the
%   simulation to the end of the pulse.
%
%   The signal, returned in S, is demodulated with the Larmor
%   frequency (w0 = gamma*B0) and is sampled at times T, specified by
%   the analog-to-digital converter object ADC.
%
%   Plotting is done at each sampling of the signal and in the
%   beginnings and ends of pulses and gradients. The best way to
%   control the plotting time step is by setting the appropriate
%   sampling interval of ADC or with the 'TimeStep' option below.
%
%   Options:
%     'Plot'       Set if plotting should be done. true (default)
%                  or false.
%     'Frame'      Set which frame to plot: 'Rotating' (default) or
%                  'Stationary'.
%     'Pause'      Pauses the simulation at times provided in an array.
%     'TimeStep'   Sets the time step of the visualization. Default
%                  determined by RF, GX, GY and (mainly) ADC. 
%     'Layout'     Select plot layout: 1 (default) or 2.
%     'PlotKSpace' Plot k-space interpretation: true or false (default).
%
%   See also IMAGINGVOLUME, RFPULSE, RECTPULSE, SINCPULSE, GRADIENT, ADC
%
%   Copyright (c) 2010-2013 Peter Nillius (nillius@mi.physics.kth.se)
