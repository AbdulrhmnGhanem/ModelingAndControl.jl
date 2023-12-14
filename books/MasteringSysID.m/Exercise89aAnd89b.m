function Exercise89a89b

% Chapter 6 Exercise 89.a and 89.b
% Exercise 89.a: Design of baseband odd and full random phase multisines with random harmonic grid
%          89.b:           bandpass
%          
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %% 
% Design of random phase multisines with a random harmonic grid %%
%                                                               %%
% Illustrated on a uniform (linear) and logarithmic frequency   %%
% distribution                                                  %%
%                                                               %%
%                                                               %%
% Rik Pintelon and Johan Schoukens                              %% 
% December 14, 2007                                             %%
%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition odd lin tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1;
fmax = 100;
fres = 2;
fs = 400;
Nblock = 3;
Spacing = 'lin';
MultiType = 'odd';
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;

% one realisation odd random harmonic grid
[ExcitedHarmOddLin, N_odd_lin, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmOddLin)

% one realisation odd random phase multisine with a given random harmonic grid
OddLinSignal = CalcMultisine(ExcitedHarmOddLin, N_odd_lin);

% plot DFT spectrum
OddLinSpec = fft(OddLinSignal)/sqrt(N_odd_lin);
SelectAll = (1:N_odd_lin/2+1).';
freqAll = (SelectAll-1)*fs/N_odd_lin;
OddLinSpec = OddLinSpec(SelectAll);
figure(1)
subplot(211)
plot(OddLinSignal);
title('Odd lin tone')
subplot(212)
plot(freqAll, db(OddLinSpec),'+')
zoom on
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition odd log tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1e2;
fmax = 1e4;
fres = 10;
% frat = 1.015;
frat = 1.1;
fs = 1.2e5;
Nblock = 3;
Spacing = 'LOG';
MultiType = 'odd';
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;
DefFreq.frat = frat;

% one realisation odd random harmonic grid
[ExcitedHarmOddLog, N_odd_log, NewDef] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmOddLog)

% one realisation odd random phase multisine with a given random harmonic grid
OddLogSignal = CalcMultisine(ExcitedHarmOddLog, N_odd_log);

% plot DFT spectrum
OddLogSpec = fft(OddLogSignal)/sqrt(N_odd_log);
SelectAll = (1:N_odd_log/2+1).';
freqAll = (SelectAll-1)*fs/N_odd_log;
OddLogSpec = OddLogSpec(SelectAll);
figure(2)
subplot(211)
plot(OddLogSignal);
title('Odd log tone')
subplot(212)
semilogx(freqAll, db(OddLogSpec),'+')
zoom on
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition full lin tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1;
fmax = 100;
fres = 1;
fs = 400;
Nblock = 3;
Spacing = 'lin';
MultiType = 'full';
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;

% one realisation odd random harmonic grid
[ExcitedHarmFullLin, N_full_lin, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmFullLin)

% one realisation odd random phase multisine with a given random harmonic grid
FullLinSignal = CalcMultisine(ExcitedHarmFullLin, N_full_lin);

% plot DFT spectrum
FullLinSpec = fft(FullLinSignal)/sqrt(N_full_lin);
SelectAll = (1:N_full_lin/2+1).';
freqAll = (SelectAll-1)*fs/N_full_lin;
FullLinSpec = FullLinSpec(SelectAll);
figure(3)
subplot(211)
plot(FullLinSignal);
title('Full lin tone')
subplot(212)
plot(freqAll, db(FullLinSpec),'+')
zoom on
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition full log tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1e2;
fmax = 1e4;
fres = 10/2;
frat = sqrt(frat);
fs = 1.2e5;
Nblock = 3;
Spacing = 'LOG';
MultiType = 'full';
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;
DefFreq.frat = frat;

% one realisation odd random harmonic grid
[ExcitedHarmFullLog, N_full_log, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmFullLog)

% one realisation odd random phase multisine with a given random harmonic grid
FullLogSignal = CalcMultisine(ExcitedHarmFullLog, N_full_log);

% plot DFT spectrum
FullLogSpec = fft(FullLogSignal)/sqrt(N_full_log);
SelectAll = (1:N_full_log/2+1).';
freqAll = (SelectAll-1)*fs/N_full_log;
FullLogSpec = FullLogSpec(SelectAll);
figure(4)
subplot(211)
plot(FullLogSignal);
title('Full log tone')
subplot(212)
semilogx(freqAll, db(FullLogSpec),'+')
zoom on
shg

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
