% Chapter 6 Exercise 95
% Prediction of the nonlinear distortions using random harmonic grid multisines
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 9 December 2010


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %%
% Prediction of the variance of the stochastic nonlinear distortions   %%
% of full and odd multisines (all harmonics are excited) via full      %%
% and odd random harmonic grid multisines.                             %% 
%                                                                      %%
% Considered system: Generalised Wiener-Hammerstein                    %%
%                                                                      %%
%       1. linear time invariant system                                %%
%       2. nonlinear FIR system                                        %%
%       3. linear time invariant system                                %%
%                                                                      %%
%                                                                      %%
% Note: no noise is added to the input/output signals                  %%
%                                                                      %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                 % sampling frequency
N = 512*2;                                              % number of data points in one period (multisine) or block (Gaussian noise)
M = 1000;                                                % number of realisations (multisine) or blocks (Gaussian noise)

% definition signal filter
cL = [1 -0.8 0.1];
dL = [1, -0.2];

% definition first LTI system
a1 = [1 -0.5 0.9];
b1 = 1;

% definition second LTI system
a2 = [1 -1.5 0.7];
b2 = [0 1 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using full multisine with all harmonics excited %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full multisine
FreqIndex_f = [1:1:N/2-1].';                                            % multisine exciting all harmonics Nyquist not included
freq_f = FreqIndex_f/N;                                                 % excited frequencies
F = length(freq_f);                                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_f);                                         % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                                     % matrix of the FRF's of all realisations

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system

% NL contribution
alpha = 0.01/20;

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex_f+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_f)));   % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex_f+1) = U(FreqIndex_f+1).*G1;                            % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = x(t) + 0.1*atan(x(t))*x(t-1)^2 + |x(t-2)|
 	z = x + alpha*atan(x) .* [x(end);x(1:end-1)].^2 + 10*alpha*abs([x(end-1:end);x(1:end-2)]);         % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                            % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex_f+1) = Z(FreqIndex_f+1).*G2;                            % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex_f+1)./U(FreqIndex_f+1)).';               % FRF of ii th realisation at the excited frequencies
end

% rms value full multisine
urms = sqrt(mean(u.^2));

G_f = (mean(Gall, 1)).';
varGs_f = (std(Gall, 0, 1).^2).';                                       % variance w.r.t. one realisation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine with all harmonics excited %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd multisine
FreqIndex_o = [1:2:N/2-1].';                                            % multisine exciting all odd harmonics Nyquist not included
freq_o = FreqIndex_o/N;                                                 % excited frequencies
F = length(freq_o);                                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_o);                                         % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                                     % matrix of the FRF's of all realisations

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex_o+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_o)));   % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                        % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex_o+1) = U(FreqIndex_o+1).*G1;                            % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = x(t) + 0.1*atan(x(t))*x(t-1)^2 + |x(t-2)|
 	z = x + alpha*atan(x) .* [x(end);x(1:end-1)].^2 + 10*alpha*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex_o+1) = Z(FreqIndex_o+1).*G2;                            % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex_o+1)./U(FreqIndex_o+1)).';               % FRF of ii th realisation at the excited frequencies
end

G_o = (mean(Gall, 1)).';
varGs_o = (std(Gall, 0, 1).^2).';                                       % variance w.r.t. one realisation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using full multisine with random harmonic grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full multisine with random harmonic grid:
% one harmonic out of Fgroup consecutive harmonics are not excited
Fgroup = 3;
FreqIndex_fr = lintone(N/2, Fgroup, 'full');
% Nyquist is not excited
if FreqIndex_fr(end) == N/2
    FreqIndex_fr(end) = [];
end % if

% excited and non-excited harmonics
MeasHarm_f = [1:N/2-1].';
NonExcitedHarm_f = HarmonicContent(MeasHarm_f, FreqIndex_fr);

freq_fr = FreqIndex_fr/N;                                               % excited frequencies
freq_frNE = NonExcitedHarm_f.inband/N;                                    % in-band non-excited frequencies
F = length(freq_fr);                                                    % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_fr);                                        % z^-1
qall = exp(-sqrt(-1)*2*pi*[0:1:N-1]/N).';
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                                     % matrix of the FRF's of all realisations
Yall_f = zeros(M, N/2-1);

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2all = polyval(fliplr(b2), qall)./polyval(fliplr(a2), qall);                    % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex_fr+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_fr)));     % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                         % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex_fr+1) = U(FreqIndex_fr+1).*G1;                                 % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = x(t) + 0.1*atan(x(t))*x(t-1)^2 + |x(t-2)|
 	z = x + alpha*atan(x) .* [x(end);x(1:end-1)].^2 + 10*alpha*abs([x(end-1:end);x(1:end-2)]);         % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                            % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y = Z.*G2all;                                                      % steady state response at all the excited frequencies
%     Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                             % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex_fr+1)./U(FreqIndex_fr+1)).';            % FRF of ii th realisation at the excited frequencies
    Yall_f(ii, :) = Y(MeasHarm_f+1).';
end

G_fr = (mean(Gall, 1)).';
varGs_fr = (std(Gall, 0, 1).^2).';                                     % variance w.r.t. one realisation

% fast method rms averaged over different realisations
varYs_f = mean(abs(Yall_f(:, NonExcitedHarm_f.inband)).^2, 1).';

% interpolation at excited frequencies
varGs_f_pred = interp1(freq_frNE, varYs_f, freq_fr)./abs(U(FreqIndex_fr+1).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine with random harmonic grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd multisine with random harmonic grid:
% one harmonic out of two consecutive odd harmonics are not excited
FreqIndex_or = lintone(N/2, Fgroup, 'odd');
% Nyquist is not excited
if FreqIndex_or(end) == N/2
    FreqIndex_or(end) = [];
end % if

% excited and non-excited harmonics
MeasHarm_o = [1:N/2-1].';     % measured harmonics
NonExcitedHarm_o = HarmonicContent(MeasHarm_o, FreqIndex_or);

freq_or = FreqIndex_or/N;                                               % excited frequencies
freq_orNE.even = NonExcitedHarm_o.even.inband/N;                        % in-band non-excited even harmonics
freq_orNE.odd = NonExcitedHarm_o.odd.inband/N;                          % in-band non-excited odd harmonics
F = length(freq_or);                                                    % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_or);                                        % z^-1
qall = exp(-sqrt(-1)*2*pi*[0:1:N-1]/N).';
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                                     % matrix of the FRF's of all realisations
Yall_o = zeros(M, N/2-1);

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2all = polyval(fliplr(b2), qall)./polyval(fliplr(a2), qall);                    % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex_or+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_or)));  % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                         % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex_or+1) = U(FreqIndex_or+1).*G1;                           % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = x(t) + 0.1*atan(x(t))*x(t-1)^2 + |x(t-2)|
 	z = x + alpha*atan(x) .* [x(end);x(1:end-1)].^2 + 10*alpha*abs([x(end-1:end);x(1:end-2)]);         % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                            % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y = Z.*G2all;                                                        % steady state response at all frequencies

    Gall(ii, :) = (Y(FreqIndex_or+1)./U(FreqIndex_or+1)).';              % FRF of ii th realisation at the excited frequencies
    Yall_o(ii, :) = Y(MeasHarm_o+1).';
end

% robust method for calculating var(Gs)
G_or = (mean(Gall, 1)).';
varGs_or = (std(Gall, 0, 1).^2).';                                     % variance w.r.t. one realisation

% fast method rms averaged over different realisations
varYs_o.even = mean(abs(Yall_o(:, NonExcitedHarm_o.even.inband)).^2, 1).'; 
varYs_o.odd = mean(abs(Yall_o(:, NonExcitedHarm_o.odd.inband)).^2, 1).'; 

% interpolation at excited frequencies
varGs_0_pred.even = interp1(freq_orNE.even, varYs_o.even, freq_or)./abs(U(FreqIndex_or+1).^2);
varGs_0_pred.odd = interp1(freq_orNE.odd, varYs_o.odd, freq_or)./abs(U(FreqIndex_or+1).^2);


%%%%%%%%%%%%%%%%%%%%%
% Compare the BLA's %
%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[50 300 1200 600])              % larger window for plotting
subplot(231)
plot(freq_f, db(G_f), 'k', freq_f, db(varGs_f)/2, 'k:', freq_o, db(G_o), 'r', freq_o, db(varGs_o)/2, 'r:')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)
title('Black: full; Red: odd', 'FontSize', 12)

subplot(232)
plot(freq_f, db(G_f), 'k', freq_f, db(varGs_f)/2, 'k:', freq_fr, db(G_fr), 'r', freq_fr, db(varGs_fr)/2, 'r:')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)
title('Black: full; Red: full-random', 'FontSize', 12)

subplot(233)
plot(freq_o, db(G_o), 'k', freq_o, db(varGs_o)/2, 'k:', freq_or, db(G_or), 'r', freq_or, db(varGs_or)/2, 'r:')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)
title('Black: odd; Red: odd-random', 'FontSize', 12)

subplot(234)
% interpolation varGs_f on odd frequencies
varGs_f_interp_o = interp1(freq_f, varGs_f, freq_o);
plot(freq_o, db(varGs_f_interp_o - varGs_o)/2, 'k')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Difference var (dB)', 'FontSize', 12)

subplot(235)
% interpolation varGs_f on full random grid
varGs_f_interp_fr = interp1(freq_f, varGs_f, freq_fr);
plot(freq_fr, db(varGs_f_interp_fr)/2 - db(varGs_fr)/2, 'k')
axis([0 0.5 0 2])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Ratio std (dB)', 'FontSize', 12)

subplot(236)
% interpolation varGs_o on odd random grid
varGs_o_interp_or = interp1(freq_o, varGs_o, freq_or);
plot(freq_or, db(varGs_o_interp_or)/2 - db(varGs_or)/2, 'k')
axis([0 0.5 0 2])
xlabel('f/fs', 'FontSize', 12)
ylabel('Ratio std (dB)', 'FontSize', 12)
       
shg
zoom on


figure(2)
set(gcf, 'Position',[50 200 1200 250])              % larger window for plotting
subplot(131)
plot(freq_fr, db(Yall_f(1, FreqIndex_fr)), 'k+', freq_frNE, db(varYs_f)/2, 'b^');
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Y(k) (dB)', 'FontSize', 12)
title('+: excited; triangle: non-excited', 'FontSize', 12)

subplot(132)
plot(freq_or, db(Yall_o(1, FreqIndex_or)), 'k+', freq_orNE.odd, db(varYs_o.odd)/2, 'ro', ...
     freq_orNE.even, db(varYs_o.even)/2, 'gx');
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Y(k) (dB)', 'FontSize', 12)
title('+: excited; x: odd non-excited; o: even non-excited', 'FontSize', 12)

subplot(133)
% interpolation varYs_o.even at odd frequencies
varYs_o_interp.even = interp1(freq_orNE.even, varYs_o.even, freq_orNE.odd);
plot(freq_fr, db(Yall_f(1, FreqIndex_fr)), 'k', freq_or, db(Yall_o(1, FreqIndex_or)), 'k--', ...
     freq_frNE, db(varYs_f)/2, 'b', freq_orNE.odd, db(varYs_o.odd + varYs_o_interp.even)/2, 'b--');
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Y(k) (dB)', 'FontSize', 12)
title('-: full; - -: odd', 'FontSize', 12)

shg;
zoom on


figure(3)
set(gcf, 'Position',[50 600 1200 250])              % larger window for plotting
subplot(131)
plot(freq_or, db(Yall_o(1, FreqIndex_or)), 'k+', freq_orNE.odd, db(varYs_o.odd)/2, 'ro', ...
     freq_orNE.even, db(varYs_o.even)/2, 'gx');
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('Y(k) (dB)', 'FontSize', 12)
title('+: excited; x: odd non-excited; o: even non-excited', 'FontSize', 12)

subplot(132)
scale = db(Fgroup/(Fgroup-1))/2;
plot(freq_f, db(G_f), 'k', freq_f, db(varGs_f)/2, 'k:', freq_o, db(varGs_o)/2, 'r:', ...
     freq_or, db(varGs_0_pred.even + varGs_0_pred.odd)/2+scale, 'k+', freq_or, db(varGs_0_pred.odd)/2+scale, 'r+')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)
title('-: robust odd; +: fast odd-random', 'FontSize', 12)

subplot(133)
plot(freq_o, db(varGs_f_interp_o - varGs_o)/2, 'k', freq_or, db(varGs_0_pred.even)/2+scale, 'k+')
axis([0 0.5 -80 40])
xlabel('f/f_s', 'FontSize', 12)
ylabel('var(G_{S, even} (dB)', 'FontSize', 12)
title('-: robust odd; +: fast odd-random', 'FontSize', 12)

shg
zoom on