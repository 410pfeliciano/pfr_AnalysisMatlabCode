function Hd = SWRFilter
%SWRFILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 10-Sep-2018 20:40:22

% FIR Window Bandpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 2000;  % Sampling Frequency

N    = 500;      % Order
Fc1  = 125;      % First Cutoff Frequency
Fc2  = 250;      % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = blackman(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
