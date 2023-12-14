% Chapter 5 Exercise 84.b
% Normality random phase multisine
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

% selection of the signal filter: see line 37
%                  phase distribution: see line 56

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 %%
% Influence of the number of frequencies on the   %%
% probability density function of a multisine     %%
%                                                 %%
%                                                 %%
% Rik Pintelon and Johan Schoukens                %% 
% April 12, 2007                                  %%
%                                                 %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                             % sampling frequency
N = 1024;                                                           % number of data points in one period (multisine) or block (Gaussian noise)
M = 10000;                                                          % number of realisations (multisine) or blocks (Gaussian noise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of the signal filter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SignalFilter = 'no';                                                % no signal filter
% SignalFilter = 'highpass';                                        % first order highpass filter
% SignalFilter = 'bandpass';                                        % 6th order bandpass filter with center frequency 0.25*fs

switch SignalFilter
    case 'no'
        cL = 1;
        dL = 1;
    case 'highpass'
        cL = [1 -0.8 0.1];
        dL = [1, -0.2];
    case 'bandpass'
        [cL, dL] = butter(3, 2*[0.20, 0.30]);
end % switch

%%%%%%%%%%%%%%%%%%%%%%%
% choice random phase %
%%%%%%%%%%%%%%%%%%%%%%%

PhaseMultisine = 'uniform';                                         % uniformly distributed in [0, pi)
PhaseMultisine = 'zeropi';                                          % zero or pi with probablity 1/2


%%%%%%%%%%%%%%%
% single sine %
%%%%%%%%%%%%%%%

steps = N/4;
FreqIndex = [steps:steps:N/2-1].';                                  % multisine exciting all harmonics Nyquist not included
freq_s = FreqIndex/N;                                               % excited frequencies
F_s = length(freq_s);                                               % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_s);                                     % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                 % filter transfer function for shaping the amplitude spectrum of the multisine

uall = zeros(N, M);
for ii = 1:M
    home, ii
    % calculate the time signal
    U = zeros(N,1);
    switch PhaseMultisine
        case 'uniform'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_s)));         % random phase multisine with prescribed amplitude spectrum
        case 'zeropi'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*pi*round(rand(size(freq_s))));	% random phase multisine with prescribed amplitude spectrum
   end
    uall(:,ii) = 2*real(ifft(U))*sqrt(N);                                           % and same rms value as the filtered Gaussian noise excitation
end % ii
stdu = std(uall(:,1));

% calculate the pdf of the multisine
[uu1, xx1_s] = hist(uall(:), sqrt(length(uall(:)))/10);               % histogram of u(t)
pdf_s = uu1/(sum(uu1)*(mean(diff(xx1_s))));                           % experimental pdf = histogram divided by its area
pdf0_s = exp(-xx1_s.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);                 % gaussian pdf with zero mean and variance stdu^2


%%%%%%%%%%%
% 3 sines %
%%%%%%%%%%%

switch PhaseMultisine
    case 'uniform'
        step = N/8;
        Number2 = '3';
    case 'zeropi'
        step = 32;
        Number2 = '15'
end
FreqIndex = [step:step:N/2-1].';                                      % multisine exciting all harmonics Nyquist not included
freq_3 = FreqIndex/N;                                                 % excited frequencies
F_3 = length(freq_3);                                                 % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_3);                                       % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                   % filter transfer function for shaping the amplitude spectrum of the multisine

uall = zeros(N, M);
for ii = 1:M
    home, ii
    % calculate the time signal
    U = zeros(N,1);
    switch PhaseMultisine
        case 'uniform'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_3)));         % random phase multisine with prescribed amplitude spectrum
        case 'zeropi'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*pi*round(rand(size(freq_3))));	% random phase multisine with prescribed amplitude spectrum
   end
    uall(:,ii) = 2*real(ifft(U))*sqrt(N);                             % and same rms value as the filtered Gaussian noise excitation
end % ii
stdu = std(uall(:,1));

% calculate the pdf of the multisine
[uu1, xx1_3] = hist(uall(:), sqrt(length(uall(:)))/10);               % histogram of u(t)
pdf_3 = uu1/(sum(uu1)*(mean(diff(xx1_3))));                           % experimental pdf = histogram divided by its area
pdf0_3 = exp(-xx1_3.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);                 % gaussian pdf with zero mean and variance stdu^2


%%%%%%%%%%%%
% 15 sines %
%%%%%%%%%%%%

switch PhaseMultisine
    case 'uniform'
        step = 2*16;
        Number3 = '15';
    case 'zeropi'
        step = 2;
        Number3 = '255'
end
FreqIndex = [step:step:N/2-1].';                                      % multisine exciting all harmonics Nyquist not included
freq_15 = FreqIndex/N;                                                % excited frequencies
F_15 = length(freq_15);                                               % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_15);                                      % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                   % filter transfer function for shaping the amplitude spectrum of the multisine
  
uall = zeros(N, M);
for ii = 1:M
    home, ii
    % calculate the time signal
    U = zeros(N,1);
    switch PhaseMultisine
        case 'uniform'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_15)));        % random phase multisine with prescribed amplitude spectrum
        case 'zeropi'
            U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*pi*round(rand(size(freq_15))));	% random phase multisine with prescribed amplitude spectrum
   end
    uall(:,ii) = 2*real(ifft(U))*sqrt(N);                                           % and same rms value as the filtered Gaussian noise excitation
end % ii
stdu = std(uall(:,1));

% calculate the pdf of the multisine
[uu1, xx1_15] = hist(uall(:), sqrt(length(uall(:)))/10);                % histogram of u(t)
pdf_15 = uu1/(sum(uu1)*(mean(diff(xx1_15))));                           % experimental pdf = histogram divided by its area
pdf0_15 = exp(-xx1_15.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);                 % gaussian pdf with zero mean and variance stdu^2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal filter for plotting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step = 2;
FreqIndex = [step:step:N/2-1].';                                        % multisine exciting all harmonics Nyquist not included
freq = FreqIndex/N;                                                     % excited frequencies
F = length(freq);                                                       % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq);                                           % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine


%%%%%%%%%%%%%%%%%%%%
% plot the figures %
%%%%%%%%%%%%%%%%%%%%

figure(1);
set(gcf, 'Position',[50 500 1200 300])                                 % larger window for plotting
subplot(141)
plot(freq, db(L))
xlabel('f/fs');
ylabel('Amplitude (dB)');
axis1 = axis;
axis([0 0.5 axis1(3:4)])
subplot(142)
plot(xx1_s, pdf_s, xx1_s, pdf0_s);
xlabel('u');
ylabel('pdf(u)');
title(['1 sine'])
subplot(143)
plot(xx1_3, pdf_3, xx1_3, pdf0_3);
xlabel('u');
ylabel('pdf(u)');
title([Number2,' sines'])
subplot(144)
plot(xx1_15, pdf_15, xx1_15, pdf0_15);
xlabel('u');
ylabel('pdf(u)');
title([Number3,' sines'])
