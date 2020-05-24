% Created with matlab R2018a by pfr 1/20/2020
% This is a bandpass filter. It uses the filtfilt function in matlab to 
% perform the filtering
% Type of filters:
%     -'theta' filter between 4-12Hz
%     -'lowGamma' filter between 25-55Hz
%     -'fastGamma' filter between 60-100Hz
%     -'ripple' filter between 125-250Hz
%     -'spindle' filter between 10-16Hz (look at Siapas and Wilson, 1998)
%     -'delta' filter between 1-4Hz
%     -'deltaSpindle' filter between 1-16Hz
%     -'highfreq' filter between 20-100Hz
% Inputs:
% -LFP is the LFP signal
% -lfpFilter has to be a string! choose between 'theta','lowGamma',
% 'fastGamma' or 'ripple'
% -Fs(Default) = 1000. Assuming that you reduced the sampling rate to 1KHz
% 
function filtLFP = bandFilterMaria (LFP, lfpFilter, Fs)

if nargin < 3
   Fs = 1000;
end
flag = 'scale';  % Sampling Flag
N    = 500;      % Order
% Create the window vector for the design algorithm.
win = blackman(N+1);

%% Filter band Selection
if strcmpi(lfpFilter, 'fastGamma') % 60-100Hz Filter
    Lowpass = 60;      
    Highpass = 100;     
elseif strcmpi(lfpFilter, 'lowGamma') % 25-55Hz Filter
    Lowpass = 25;     
    Highpass = 55;      
elseif strcmpi(lfpFilter, 'theta') % 4-12Hz Filter
    Lowpass = 4;      
    Highpass = 12;    
elseif strcmpi(lfpFilter, 'ripple') % 100-275Hz Filter
    Lowpass = 125;
    Highpass = 250;
elseif strcmpi(lfpFilter, 'spindle') % 10-16Hz Filter
    Lowpass = 10;
    Highpass = 16;
elseif strcmpi(lfpFilter, 'delta') % 1-4Hz Filter
    Lowpass = 1;
    Highpass = 4;
elseif strcmpi(lfpFilter, 'deltaSpindle') % 1-16Hz Filter
    Lowpass = 1;
    Highpass = 16;
elseif strcmpi(lfpFilter, 'highFreq') % 100-450Hz Filter
    Lowpass = 100;
    Highpass = 375;
elseif strcmpi(lfpFilter, 'sharpwave') % 1-4Hz Filter
    Lowpass = 1;
    Highpass = 50;
end
%%

% Calculate the coefficients using the FIR1 function.
fprintf(1, 'Filtering the data\n');
Rb  = fir1(N, [Lowpass Highpass]/(Fs/2), 'bandpass', win, flag);
        bandFiltHd = dfilt.dffir(Rb); % Ripple Filter
        filtLFP = filtfilt(bandFiltHd.Numerator,1,LFP);
 fprintf(1, 'Done!\n'); 
 %%
%   prompt = 'Assign a new variable name to your data ? ';
%   xx = input(prompt, 's');
%   assignin('base',xx,filtLFP)
% %% saving
% if nargin > 3
%     prompt = 'What is the channel # ? ';
%     channeltosave = input(prompt);
%     fprintf(1, 'Saving Files!\n');
%     file_name1 = sprintf('ch%dLFP%s.mat',channeltosave,lfpFilter);
%     save (file_name1, xx,'-v7.3') % saving spike data
%     fprintf(1, 'Done!\n'); 
% end
end