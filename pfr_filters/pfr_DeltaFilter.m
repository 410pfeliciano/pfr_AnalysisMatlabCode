function deltaLFP = pfr_DeltaFilter(LFP, Fs, varargin)
%{
DELTAFILTER Returns a discrete-time filter object.
Created with matlab R2018a BY pfr 1/20/2020
FIR Window Lowpass filter designed using the FIR1 function.
Inputs
1- LFP - lfp data
2-Fs(default) = 1000 if you downsampled the data to 1KHz
3-VARARGIN - Put any number/string if you want to save the data into the 
working directory
%}
%All frequency values are in Hz.
if nargin < 2
   Fs = 1000;
end
N    = 500;      % Order
Fc   = 4;        % Cutoff Frequency
flag = 'scale';  % Sampling Flag

% Create the window vector for the design algorithm.
win = blackman(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);
fprintf(1, 'Filtering the data\n');
deltaLFP = filtfilt(Hd.Numerator,1,LFP);
fprintf(1, 'Filtering Done!\n'); 
%%
if nargin > 2
    fprintf(1, 'Saving Files!\n');
    file_name1 = sprintf('deltaLFP.mat');
    save (file_name1, 'deltaLFP','-v7.3') % saving spike data
    fprintf(1, 'Done!\n'); 
end
end