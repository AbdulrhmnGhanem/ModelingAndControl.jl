% Chapter 5 Exercise 84.c
% Influence of harmonic content multisine on BLA
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %%
% Comparison BLA obtained by full, odd, odd-odd, odd-random   %%
% (harmonic grid), odd-sparse (harmonic grid) multisines with %%
% the same power spectrum (colouring and rms value)           %%
%                                                             %%
% Considered system: Generalised Wiener-Hammerstein           %%
%                                                             %%
%       1. linear time invariant system                       %%
%       2. nonlinear FIR system                               %%
%       3. linear time invariant system                       %%
%                                                             %%
%                                                             %%
% Rik Pintelon and Johan Schoukens                            %% 
% April 12, 2007                                              %%
% version May 21, 2007                                        %%
%                                                             %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                 % sampling frequency
N = 512*2;                                              % number of data points in one period (multisine) or block (Gaussian noise)
M = 10000;                                              % number of realisations (multisine) or blocks (Gaussian noise)

% definition signal filter
cL = [1 -0.8 0.1];
dL = [1, -0.2];

% definition first LTI system
a1 = [1 -0.5 0.9];
b1 = 1;

% definition second LTI system
a2 = [1 -1.5 0.7];
b2 = [0 1 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using full multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full multisine
FreqIndex = [1:1:N/2-1].';                                              % multisine exciting all harmonics Nyquist not included
freq_f = FreqIndex/N;                                                   % excited frequencies
F = length(freq_f);                                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_f);                                         % z^-1
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_f)));     % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

% rms value full multisine
urms = sqrt(mean(u.^2));

G_f = (mean(Gall, 1)).';
varG_f = (std(Gall, 0, 1).^2/M).';                                      % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd multisine
FreqIndex = [1:2:N/2-1].';                                              % multisine exciting all odd harmonics Nyquist not included
freq_o = FreqIndex/N;                                                   % excited frequencies
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_o)));     % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                        % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_o = (mean(Gall, 1)).';
varG_o = (std(Gall, 0, 1).^2/M).';                                      % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd-odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd-odd multisine
FreqIndex = [1:4:N/2-1].';                                              % multisine exciting all odd harmonics Nyquist not included
freq_oo = FreqIndex/N;                                                  % excited frequencies
F = length(freq_oo);                                                    % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_oo);                                        % z^-1
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_oo)));     % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                                  % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_oo = (mean(Gall, 1)).';
varG_oo = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine with random harmonic grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd multisine with random harmonic grid:
% one harmonic out of two consecutive odd harmonics are not excited
FreqIndex = lintone(N/2, 2);
freq_or = FreqIndex/N;                                                  % excited frequencies
F = length(freq_or);                                                    % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_or);                                        % z^-1
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_or)));     % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                                  % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_or = (mean(Gall, 1)).';
varG_or = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using lacunar odd multisines %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd-sparse multisine: 16 frequencies
FreqIndex = [1:32:N/2-1].';                                             % multisine exciting all odd harmonics Nyquist not included
freq_o_16 = FreqIndex/N;                                                % excited frequencies
F = length(freq_o_16);                                                  % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_o_16);                                      % z^-1
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_o_16)));  % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                        % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_o_16 = (mean(Gall, 1)).';
varG_o_16 = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


% lacunar odd multisine: 32 frequencies
FreqIndex = [1:16:N/2-1].';                                             % multisine exciting all odd harmonics Nyquist not included
freq_o_32 = FreqIndex/N;                                                % excited frequencies
F = length(freq_o_32);                                                  % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_o_32);                                      % z^-1
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
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_o_32)));  % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                        % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_o_32 = (mean(Gall, 1)).';
varG_o_32 = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


% lacunar odd multisine: 64 frequencies
FreqIndex = [1:8:N/2-1].';                                             % multisine exciting all odd harmonics Nyquist not included
freq_o_64 = FreqIndex/N;                                               % excited frequencies
F = length(freq_o_64);                                                 % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_o_64);                                     % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                    % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                                    % matrix of the FRF's of all realisations

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                   % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                   % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_o_64))); % random phase multisine with prescribed amplitude spectrum
    u = real(ifft(U));
    u = u/sqrt(mean(u.^2))*urms;                                       % and same rms value as the full multisine excitation
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                               % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2 + 10*|x(t-2)|
 	z = atan(x) .* [x(end);x(1:end-1)].^2 + 10*abs([x(end-1:end);x(1:end-2)]);      % since x is periodic x(t-1) and x(t-2) are obtained
                                                                                    % via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_o_64 = (mean(Gall, 1)).';
varG_o_64 = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


%%%%%%%%%%%%%%%%%%%%%
% Compare the BLA's %
%%%%%%%%%%%%%%%%%%%%%

% comparison full and odd multisines
figure(1)
set(gcf, 'Position',[50 40 1200 900])              % larger window for plotting
subplot(331)
plot(freq_f, db(G_f), 'k', freq_f, db(varG_f)/2, 'k:', freq_o, db(G_o), 'r', freq_o, db(varG_o)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(332)
plot(freq_f, unwrap(angle(G_f))*180/pi, 'k', freq_o, unwrap(angle(G_o))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: full multisine; Red: odd multisine', 'FontSize', 12)

subplot(333)
plot(freq_o, db(G_f(1:2:end)-G_o), 'k', freq_o, db(varG_f(1:2:end) + varG_o)/2, 'r');                  % variance(Gf-Go) = var(Gf)+var(Go)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12)

% comparison odd and odd-odd multisines
subplot(334)
plot(freq_o, db(G_o), 'k', freq_o, db(varG_o)/2, 'k:', freq_oo, db(G_oo), 'r', freq_oo, db(varG_oo)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(335)
plot(freq_o, unwrap(angle(G_o))*180/pi, 'k', freq_oo, unwrap(angle(G_oo))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: odd multisine; Red: odd-odd multisine', 'FontSize', 12)

subplot(336)
plot(freq_oo, db(G_o(1:2:end)-G_oo), 'k', freq_oo, db(varG_o(1:2:end) + varG_oo)/2, 'r');                  % variance(Go-Goo) = var(Go)+var(Goo)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12)
       
% comparison odd and odd-random multisines
subplot(337)
plot(freq_o, db(G_o), 'k', freq_o, db(varG_o)/2, 'k:', freq_or, db(G_or), 'r', freq_oo, db(varG_or)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(338)
plot(freq_o, unwrap(angle(G_o))*180/pi, 'k', freq_or, unwrap(angle(G_or))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: odd multisine; Red: odd-random multisine', 'FontSize', 12)

subplot(339)
% interpolation odd multisine measurements to have it at the same
% frequency grid as the odd-random multisine measurements
G_o_interp = interp1(freq_o, G_o, freq_or);
varG_o_interp = interp1(freq_o, varG_o, freq_or);
plot(freq_or, db(G_o_interp - G_or), 'k', freq_or, db(varG_o_interp + varG_or)/2, 'r');                  % variance(Go-Goo) = var(Go)+var(Goo)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12)

zoom on;
shg

% comparison variance odd and odd-odd multisines
figure(2);
plot(freq_o, db(varG_f(1:2:end))/2 - db(varG_o)/2, 'k', freq_oo, db(varG_o(1:2:end))/2 - db(varG_oo)/2, 'r',...
     freq_or, db(varG_o_interp)/2 - db(varG_or)/2, 'r:');
axis([0 0.5 0 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('\Delta std(G_{BLA}) (dB)', 'FontSize', 12)
title('Black: full - odd; Red: diff. odd - odd-odd; Red --: diff. odd - odd-random', 'Fontsize', 12)

% comparison odd and lacunar-odd multisines
figure(3)
set(gcf, 'Position',[50 40 1200 900])              % larger window for plotting

% comparison odd multisine and odd multisine with 16 freq.
Select = [1:16:length(freq_o)];
subplot(331)
plot(freq_o(Select), db(G_o(Select)), 'k', freq_o(Select), db(varG_o(Select))/2, 'k:', freq_o_16, db(G_o_16), 'r', freq_o_16, db(varG_o_16)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(332)
plot(freq_o(Select), unwrap(angle(G_o(Select)))*180/pi, 'k', freq_o_16, unwrap(angle(G_o_16))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: 256 frequencies, and gray: 16 frequencies', 'FontSize', 12)

subplot(333)
plot(freq_o_16, db(G_o(Select)-G_o_16), 'k', freq_o_16, db(varG_o(Select) + varG_o_16)/2, 'r');                  % variance(Gf-Go) = var(Gf)+var(Go)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12)

% comparison odd multisine and odd multisine with 32 freq.
Select = [1:8:length(freq_o)];
subplot(334)
plot(freq_o(Select), db(G_o(Select)), 'k', freq_o(Select), db(varG_o(Select))/2, 'k:', freq_o_32, db(G_o_32), 'r', freq_o_32, db(varG_o_32)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(335)
plot(freq_o(Select), unwrap(angle(G_o(Select)))*180/pi, 'k', freq_o_32, unwrap(angle(G_o_32))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: 256 frequencies, and gray: 32 frequencies', 'FontSize', 12)

subplot(336)
plot(freq_o_32, db(G_o(Select)-G_o_32), 'k', freq_o_32, db(varG_o(Select) + varG_o_32)/2, 'r');                  % variance(Go-Go_32) = var(Go)+var(Go_32)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12)

% comparison odd multisine and odd multisine with 64 freq.
Select = [1:4:length(freq_o)];
subplot(337)
plot(freq_o(Select), db(G_o(Select)), 'k', freq_o(Select), db(varG_o(Select))/2, 'k:', freq_o_64, db(G_o_64), 'r', freq_o_64, db(varG_o_64)/2, 'r:')
axis([0 0.5 -50 50])
xlabel('f/fs', 'FontSize', 12)
ylabel('G_{BLA} (dB)', 'FontSize', 12)

subplot(338)
plot(freq_o(Select), unwrap(angle(G_o(Select)))*180/pi, 'k', freq_o_64, unwrap(angle(G_o_64))*180/pi, 'r')
axis([0 0.5 -400 0])
xlabel('f/fs', 'FontSize', 12)
ylabel('arg(G_{BLA}) (o)', 'FontSize', 12)
title('Black: 256 frequencies, and gray: 64 frequencies', 'FontSize', 12)

subplot(339)
plot(freq_o_64, db(G_o(Select)-G_o_64), 'k', freq_o_64, db(varG_o(Select) + varG_o_64)/2, 'r');                  % variance(Go-Go_64) = var(Go)+var(Go_64)
axis([0 0.5 -60 20])
xlabel('f/fs', 'FontSize', 12, 'FontSize', 12)
ylabel('error (dB)', 'FontSize', 12, 'FontSize', 12)

shg
zoom on

