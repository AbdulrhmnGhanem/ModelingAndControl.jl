function Exercise89c

% Chapter 6 Exercise 89.c
% Fast method for noisy input/output measurements - open loop example
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 8 December 2010

%  choice of the type of multisine (full, odd, random grid)
%       at lines: 52-56
%  lin or log resolution: at line 42-47

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %% 
% Simulation of a Wiener-Hammerstein system, consisting of      %%
% cascade of a second order discrete-time system, a polynomial  %%
% static nonlinearity, and a second order discrete-time system. %%
%                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition excitation signal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 700e6;
fmax = 1100e6;
fres = 1e6;
fs = 4e9;
frat = 1.02;
Nblock = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choice frequency spacing %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Spacing = 'lin';
%    Spacing = 'log'; fmin = 10e6;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choice type multisine  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    TypeMultisine = 'odd';
    OddOdd = true;
    OddOdd = false;
%     TypeMultisine = 'full';



if strcmp(TypeMultisine, 'full')
    fres = fres/2;
end % if


DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fs = fs;
DefFreq.fres = fres;
DefFreq.frat = frat;
[ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMultisine);

% odd-odd multisine
if OddOdd
    StarHar = round(fmin/(fs/N));
    StopHar = round(fmax/(fs/N));
    ExcitedHarm = ([StarHar+1:4:StopHar-1]).';
end % if odd-odd


%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation measurement %
%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 20;                     % number of realisations random phase multisine
% M = 1;
P = 6;                      % number of consecutive periods

freqall = ([0:1:N-1]*fs/N).';               % all DFT frequencies over full unit circle
qall = exp(-sqrt(-1)*2*pi*freqall/fs);      % z^-1 over full unit circle
qOdE = exp(-sqrt(-1)*2*pi*ExcitedHarm/N);   % z^-1 at odd excited frequencies
fstart = 600e6;                             % first measured frequency
fstart = 1e6;
if Spacing == 'log'
    if strcmp(TypeMultisine, 'odd')
        fstart = fres / 2;
    else % full multisine
        fstart = fres ;
    end % if
    fstop = fmax * 3;
end % if
fstop = 1200e6;                             % last measured frequency
fstop = 2000e6;                         % last measured frequency
if strcmp(TypeMultisine, 'odd')
    freqmeas = ([fstart:fres/2:fstop]).';       % all measured frequencies
else % full multisine
    freqmeas = ([fstart:fres:fstop]).';       % all measured frequencies
end % if
MeasHarm = freqmeas/(fs/N);                 % measured harmonics
Fall = length(MeasHarm);                    % number of measured frequencies
Uall = zeros(M, P, Fall);                   % all input spectra at the measured frequencies
Yall = zeros(M, P, Fall);                   % all output spectra at the measured frequencies


% discrete-time LTI1
a1 = [1, -0.2, 0.9];
b1 = [1];
G1all = polyval(fliplr(b1), qall)./polyval(fliplr(a1), qall);
G1 = polyval(fliplr(b1), qOdE)./polyval(fliplr(a1), qOdE);

% polynomial static nonlinearity
alfa = 0.1;                                % coefficient x^2
beta = 0.001;                              % coefficient x^3
if Spacing == 'log'
    beta = 0.1;
end % if

% discrete-time LTI2
a2 = [1, -0.5, 0.9];
b2 = [1, 0.5];
G2all = polyval(fliplr(b2), qall)./polyval(fliplr(a2), qall);
G2 = polyval(fliplr(b2), qOdE)./polyval(fliplr(a2), qOdE);

% input/output noise standard deviations
sigmau = 0.1;      % input noise std
sigmay = 0.2;      % output noise std

for ii=1:M
    
    % generation one period random phase multisine with rms value one
    u = CalcMultisine(ExcitedHarm, N);
    
   % response to the first LTI system
    U = fft(u)/sqrt(N);
    X = U.*G1all;
    x = real(ifft(X)*sqrt(N));

    % response of the static nonlinearity
    z = x + alfa*x.^2 + beta*x.^3;

    % response to the second LTI system
    Z = fft(z)/sqrt(N);
    Y = Z.*G2all;
    
    % collect measurements
    Uall(ii,:,:) = repmat((U(MeasHarm+1)).', P, 1);
    Yall(ii,:,:) = repmat((Y(MeasHarm+1)).', P, 1);
   
end

% noiseless input
U0all = Uall;
% add white input and output noise to the measurements
Uall = Uall + (randn(size(Uall)) + sqrt(-1)*randn(size(Uall)))* sigmau/sqrt(2);
Yall = Yall + (randn(size(Yall)) + sqrt(-1)*randn(size(Yall))) * sigmay/sqrt(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required variables:                                                    %
%   ExcitedHarm     =   excited odd harmonics in multisine               %
%                           size: F x 1                                  %
%   MeasHarm        =   all measured harmonics                           %
%                           size: Fall x 1                               %
%   N               =   number of samples in one period of the multisine %
%   fs              =   sampling frequency                               %
%   Uall, Yall      =   M x P X Fall input output spectra where          %
%                           M = number of realisations multisine         %
%                           P = number of consecutive periods multisine  %
%                           Fall = number of all measured frequencies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of all measured harmonics
MeasHarm = round(freqmeas/(fs/N));

% NL analysis on IO spectra
[Y, Yc, U, G, freq] = Fast_NL_Anal(Yall, Uall, ExcitedHarm, MeasHarm, fs, N);

% plot results
RealNum = 5;            % number of realisation to be plotted
FigNum = 1;             % number of first figure
Plot_Fast_NL_Anal(Y, Yc, U, G, freq, Spacing, RealNum, FigNum);



function [ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMulti);
%
%       function [ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMulti);
%
%       calculates the harmonic content of an odd or full random harmonic grid multisine signal
%       with a linear or a logarithmic spacing in the band defined by [fmin, fmax]
%       for a lowpass design
%           fmin = fres/2	for odd multisines
%           fmin = fres     for full multisines
%       where
%           fres =          the frequency spacing between two consecutive odd harmonics for odd multisines 
%           fres =          the frequency spacing between two consecutive harmonics for full multisines 
%
%
%   OUTPUT
%
%       ExcitedHarm     =   excited frequencies of the multisine expressed in harmonic numbers
%
%       N               =   number of time domain samples in one period
%
%       NewDefFreq      =   struct containing the information about the sampling frequency, frequency resolution,
%                           and the excited frequency band {'fs', 'fmin', 'fmax', 'fres', 'frat'}
%                               NewDefFreq.fs      =    sampling frequency of the generator (remains unchanged!)
%                               NewDefFreq.fres    =    frequency spacing in Hz between the (odd) harmonics 
%                                                       fres is modified such that fs/fres is an integer number
%                                                       Note: not all (odd) harmonics are excited; especially for the log spacing
%                               NewDefFreq.fmin    =    lowest excited frequency in Hz
%                                                       fmin is modified such that fmin/fres is an (odd) integer number
%                                                       (for lowpass design fmin = fres/2)
%                               NewDefFreq.fmax    =    largest excited frequency in Hz
%                                                       fmax is modified such that fmax/fres is an (odd) integer number
%
%   INPUT
%
%       DefFreq         =   struct containing the information about the sampling frequency, frequency resolution,
%                           and the excited frequency band {'fs', 'fmin', 'fmax', 'fres', 'frat'}
%                               DefFreq.fs      =   sampling frequency of the generator
%                               DefFreq.fres    =   frequency spacing in Hz between the (odd) harmonics 
%                               DefFreq.fmin    =   lowest excited frequency in Hz
%                                                   (for lowpass design of odd multisines fmin = fres/2)
%                               DefFreq.fmax    =   largest excited frequency in Hz
%                               DefFreq.frat    =   ratio between consecutive (odd) harmonics for a logarithmic frequency spacing
%
%       Nblock          =   size of the group of consecutive (odd) harmonics where one is randomly eliminated
%                           if Nblock = Inf then no detection lines are present 
%
%       Spacing         =   linear or logarithmic frequency spacing; optional parameter (default = 'linear') 
%                               'lin', 'linear':       linear frequency spacing
%                               'log', 'logarithmic':  quasi logarithmic frequency spacing (rounded to DFT grid)
%
%       TypeMulti       =   type multisine; optional (default = 'odd')
%                               'odd':   odd random harmonic grid multisine        
%                               'full':  full random harmonic grid multisine 
%
% Rik Pintelon, March 2006
% version December 4, 2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = DefFreq.fs;
fres = DefFreq.fres;
fmin = DefFreq.fmin;
fmax = DefFreq.fmax;

if isfield(DefFreq, 'frat')
	if ~isempty(DefFreq.frat)
        frat = DefFreq.frat;
	end % if not empty
end % if isfield

% type of the frequency spacing
if nargin == 2
    Spacing = 'lin';
end % if
Spacing = lower(Spacing);

% type of the multisine
if nargin <= 3
    TypeMulti = 'odd';
end % if
TypeMulti = lower(TypeMulti);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation odd random harmonic grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(TypeMulti, 'odd')
    fres = fres/2;                  % frequency spacing between odd and even harmonics
end % if odd multisine
N = ceil(fs/fres);                  % number of time domain points in one period
fres = fs/N;
FreqSpan = floor(fmax/fres);        % frequency span for a lowpass signal

% calculate the (odd) random harmonic numbers for a lowpass signal
switch Spacing
    
    case {'lin', 'linear'}				
		ExcitedHarm = lintone(FreqSpan, Nblock, TypeMulti);
		
    case {'log', 'logarithmic'}
        ExcitedHarm = logtone(FreqSpan, frat, Nblock, TypeMulti);
        
end % switch Spacing

% convert lowpass in bandpass
FirstNonZeroHarm = ceil(fmin/fres);
RemoveHarm = find(ExcitedHarm < FirstNonZeroHarm);
ExcitedHarm(RemoveHarm) = [];


%%%%%%%%%%%%%%%%%%%%%%%%
% realised frequencies %
%%%%%%%%%%%%%%%%%%%%%%%%

NewDefFreq.fs = fs;
NewDefFreq.fres = fres;
if strcmp(TypeMulti, 'odd')
    NewDefFreq.fres = 2*NewDefFreq.fres;       % frequency spacing between the odd harmonics
end % if odd multisine
NewDefFreq.fmin = ExcitedHarm(1)*fres;
NewDefFreq.fmax = ExcitedHarm(end)*fres;

