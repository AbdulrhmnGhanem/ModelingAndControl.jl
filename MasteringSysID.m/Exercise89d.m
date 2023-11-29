function Exercise89d

% Chapter 6 Exercise 89.d
% Fast method for noisy input/output measurements - closed loop example
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 8 December 2010

%  choice of the type of multisine: lines 45 - 52


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %% 
% Simulation of a Wiener system within a unit delay feedback. %%
% The Wiener system consists of a first order discrete-time   %%
% system followed by the static nonlinearity 'tanh'           %%
%                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition excitation signal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 0.01;
fmax = 2;
fres = 0.01;
fs = 160;
Nblock = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choice frequency spacing %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Spacing = 'lin';

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

[ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMultisine);

% odd-odd multisine
if OddOdd
    StarHar = round(fmin/(fs/N));
    StopHar = round(fmax/(fs/N));
    ExcitedHarm = ([StarHar+1:4:StopHar-1]).';
end % if odd-odd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition Wiener system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialisation of the structures
a = struct('one', [], 'two', []);
b = a;

% first discrete-time
tau = 0.5;
a.one = [1, -exp(-1/(tau*fs))];
b.one = 100*(1-exp(-1/(tau*fs)));

% second discrete-time
a.two = 1;
b.two = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation measurement %
%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 20;                     % number of realisations random phase multisine
% M = 1;
P = 6;                      % number of consecutive periods

fstart = fmin;                             % first measured frequency
fstop = 2*fmax;
if strcmp(TypeMultisine, 'odd')
    freqmeas = ([fstart:fres/2:fstop]).';       % all measured frequencies
else % full multisine
    freqmeas = ([fstart:fres:fstop]).';       % all measured frequencies
end % if
MeasHarm = round(freqmeas/(fs/N));                 % measured harmonics
Fall = length(MeasHarm);                    % number of measured frequencies
Uall = zeros(M, P, Fall);                   % all input spectra at the measured frequencies
Yall = zeros(M, P, Fall);                   % all output spectra at the measured frequencies

% initialisation of the iteration variables
% for solving the steady state response of the
% unit delay feedback loop
ItterVar = struct('MaxNumber', [], 'RelError', []);
ItterVar.MaxNumber = 100;
ItterVar.RelError = 1e-12;

% input/output noise standard deviations
sigmau = 0.1/1000;      % input noise std
sigmay = 0.2/1000;      % output noise std

% rms value input
alpha = 0.2;

for ii=1:M
    
    ii
    % generation one period random phase multisine with rms value alpha
    r = alpha*CalcMultisine(ExcitedHarm, N);
    
    % one sample delay feedback
    [y0, u0] = WHfeedback(b, a, 'tanh', r, ItterVar);

    % input/output DFT spectra
    U0  = fft(u0)/sqrt(N);
    Y0  = fft(y0)/sqrt(N);
    
    % collect measurements
    Uall(ii,:,:) = repmat((U0(MeasHarm+1)).', P, 1);
    Yall(ii,:,:) = repmat((Y0(MeasHarm+1)).', P, 1);
   
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
[Y0, Yc0, U0] = Fast_NL_Anal(Y0(MeasHarm+1).', U0(MeasHarm+1).', ExcitedHarm, MeasHarm, fs, N);

% plot results
RealNum = 5;            % number of realisation to be plotted
FigNum = 1;             % number of first figure
Plot_Fast_NL_Anal(Y, Yc, U, G, freq, Spacing, RealNum, FigNum);

% calculation -1/feedback transfer function
figure;
switch TypeMultisine
    case 'odd'
        minusinvGfeeback = mean(Y0.mean.NE.odd.all,1)./mean(U0.mean.NE.odd.all, 1);
        freqNEall = freq.NE.odd.all;
    case 'full'
        minusinvGfeeback = mean(Y0.mean.NE.all,1)./mean(U0.mean.NE.all, 1);        
        freqNEall = freq.NE.all;
end % switch
% true -1/feedback transfer function
H0 = -exp(sqrt(-1)*2*pi*freqNEall/fs).';
subplot(211)
plot(freqNEall, db(minusinvGfeeback), freqNEall, db(minusinvGfeeback-H0))
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
title('FRF at the non-excited harmonics')
subplot(212)
plot(freqNEall, unwrap(angle((minusinvGfeeback)))*180/pi)
xlabel('Frequency (Hz)')
ylabel('Phase (o)')
title('FRF at the non-excited harmonics')



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



function [y, u, RelError] = WHfeedback(b, a, fNL, r, Itter);
%
%       function [y, u, RelError] = WHfeedback(b, a, fNL, r, Itter);
%
%       Calculates the response of a discrete-time Wiener-Hammerstein
%       system within a one sample delay feedback to a periodic input. 
%
%
%   Output parameters
%
%       y           =   output Wiener-Hammerstein system 
%       u           =   input Wiener-Hammerstein system
%       RelError    =   actual relative error of the calculations
%
%
%   Input parameters
%
%       b           =   structure ('one', 'two') containing the numerator coefficients of the linear dynamic systems
%                       in increasing powers of z^-1
%                           b.one	=   numerator coefficients first linear dynamic system 
%                           b.two	=   numerator coefficients second linear dynamic system 
%
%       a           =   structure ('one', 'two') containing the denominator coefficients of the linear dynamic systems
%                       in increasing powers of z^-1
%                           a.one	=   denominator coefficients first linear dynamic system 
%                           a.two	=   denominator coefficients second linear dynamic system 
%
%       fNL         =   string containing the (m-file) function describing the static nonlinear block 
%
%       r           =   one period of the periodic input of the feedback system 
%
%       Itter       =   structure {'MaxNumber', 'RelError'} containing the itteration parameters. 
%                           Itter.MaxNumber     =   maximal number of itterations 
%                           Itter.RelError      =   maximal relative error of the calculations 
%
%   Rik Pintelon
%   January 9, 2008
%

% normalisation leading denominator coefficients of the linear dynamic systems
b.one = b.one(:).';
a.one = a.one(:).';
b.two = b.two(:).';
a.two = a.two(:).';
b.one = b.one/a.one(1);
a.one = a.one/a.one(1);
b.two = b.two/a.two(1);
a.two = a.two/a.two(1);

% initialisation of the variables
nb1 = length(b.one)-1;
nb2 = length(b.two)-1;
na1 = length(a.one)-1;
na2 = length(a.two)-1;
r = r(:);
N = length(r);
e = zeros(N+nb1, 1);
w = zeros(N+na1, 1);
z = zeros(N+nb2, 1);
y = zeros(N+na2+1, 1);
RelError = inf;
NumberOfItter = 0;

while (RelError > Itter.RelError) & (NumberOfItter <= Itter.MaxNumber)
    
    yold = y;
    NumberOfItter = NumberOfItter + 1;
    e(1:nb1) = e(end-nb1+1:end);
    w(1:na1) = w(end-na1+1:end);
    z(1:nb2) = z(end-nb2+1:end);
    y(1:na2+1) = y(end-na2:end);

    % recursive calculation of the output
    for tt = 1:N

        % error signal
        e(tt+nb1) = r(tt) - y(tt+na2);

        % output first dynamic system
        w(tt+na1) = b.one * flipud(e(tt:tt+nb1)) - a.one(2:end) * flipud(w(tt:tt+na1-1));

        % output static nonlinearity
        z(tt+nb2) = feval(fNL, w(tt+na1));

        % output second dynamic system
        y(tt+na2+1) = b.two * flipud(z(tt:tt+nb2)) - a.two(2:end) * flipud(y(tt+1:tt+na2));

    end % tt
    
    % relative variation output
    RelError = max(abs((y(na2+2:end) - yold(na2+2:end))./y(na2+2:end)));

end % while

% remove the first na2+1 output samples
y = y(na2+2:end);
u = e(nb1+1:end);
