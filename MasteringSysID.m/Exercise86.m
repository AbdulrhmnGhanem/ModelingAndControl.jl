% Chapter 5 Exercise 86
% BLA of a cascade
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %%
% Comparison BLA of the cascade of two nonlinear systems  %%
% with the cascade of the BLA's of each system seperately %%
%                                                         %%
%                                                         %%
% Rik Pintelon and Johan Schoukens                        %% 
% April 16, 2007                                          %%
%                                                         %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                 % sampling frequency
N = 512*2;                                              % number of data points in one period (multisine) or block (Gaussian noise)
M = 40000;                                              % number of realisations (multisine) or blocks (Gaussian noise)

% definition signal filter
cL = [1 -0.8 0.1];
dL = [1, -0.2];

% definition first LTI system
a1 = [1 -0.5 0.9];
b1 = 1;

% definition second LTI system[b1, a1] = cheby1(7, 6, 2*0.3);
a2 = [1 -1.5 0.7];
b2 = [0 1 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full multisine
FreqIndex = [1:2:N/2-1].';                                              % multisine exciting all odd harmonics Nyquist not included
freq = FreqIndex/N;                                                     % excited frequencies
F = length(freq);                                                       % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq);                                           % z^(-1)
qf = exp(-sqrt(-1)*2*pi*[0:N-1]/N).';                                   % full unit circle
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);                     % filter transfer function for shaping the amplitude spectrum of the multisine
G1all = zeros(M, F);                                                    % matrix FRF's of all realisations of the first NL system
G2all = zeros(M, F);                                                    % matrix FRF's of all realisations of the second NL system
Gall = zeros(M, F);                                                     % matrix FRF's of the cascade of the two NL systems

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system
G1f = polyval(fliplr(b1), qf)./polyval(fliplr(a1), qf);

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system
G2f = polyval(fliplr(b2), qf)./polyval(fliplr(a2), qf);

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq)));       % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);
    U = fft(u)/sqrt(N);
    
    %%%%%%%%%%%%%%%%%%%%
    % First GWH system %
    %%%%%%%%%%%%%%%%%%%%

    % first LTI system
    X1 = zeros(N,1);
	X1(FreqIndex+1) = U(FreqIndex+1).*G1;                               % steady state response at the excited frequencies
    x1 = 2*real(ifft(X1))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2
 	z1 = x1 + 0.5*atan(x1) .* [x1(end);x1(1:end-1)].^2;                 % since x1 is periodic x1(t-1) is obtained
                                                                        % via a circular shift
    Z1 = fft(z1)/sqrt(N);
    
	% second LTI system
    Y1 = zeros(N,1);
    Y1 = Z1.*G2f;                                                       % steady state response at all frequencies

    G1all(ii, :) = (Y1(FreqIndex+1)./U(FreqIndex+1)).';                 % FRF of ii th realisation at the excited frequencies
    
    
    %%%%%%%%%%%%%%%%%%%%
    % Second WH system %
    %%%%%%%%%%%%%%%%%%%%

    % first LTI system
    X2 = zeros(N,1);
	X2 = Y1.*G1f;                                                       % steady state response at all frequencies
    x2 = real(ifft(X2))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t))*x(t-1)^2
 	z2 = x2 + 0.5*atan(x2) .* [x2(end);x2(1:end-1)].^2;                 % since x2 is periodic x2(t-1) is obtained
    Z2 = fft(z2)/sqrt(N);
    
	% second LTI system
    Y2 = zeros(N,1);
    Y2 = Z2.*G2f;                                                       % steady state response at all frequencies

    G2all(ii, :) = (Y2(FreqIndex+1)./Y1(FreqIndex+1)).';                % FRF of ii th realisation at the excited frequencies
    
    
    %%%%%%%%%%%%%%%
    % FRF cascade %
    %%%%%%%%%%%%%%%
    
    Gall(ii, :) = (Y2(FreqIndex+1)./U(FreqIndex+1)).';                  % FRF of ii th realisation at the excited frequencies
    
end

% BLA first NL system
Gbla1 = (mean(G1all, 1)).';
varGbla1 = (std(G1all, 0, 1).^2/M).';                                   % variance mean value

% BLA second NL system
Gbla2 = (mean(G2all, 1)).';
varGbla2 = (std(G2all, 0, 1).^2/M).';                                   % variance mean value

% variance cascade BLA's (the difficulty is that Ys2 is not
% uncorrelated with Ys1)
G1all = (G1all - repmat(Gbla1.', M, 1)) .* repmat(Gbla2.', M, 1);
G2all = (G2all - repmat(Gbla2.', M, 1)) .* repmat(Gbla1.', M, 1);
G1all = G1all + G2all;
varG1G2 = (std(G1all, 0, 1).^2/M).';

% BLA cascade systen
Gbla = (mean(Gall, 1)).';
varGbla = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison BLA(cascade) with cascade(BLA) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[50 500 1200 300])              % larger window for plotting
subplot(131)
plot(freq, db(Gbla1.*Gbla2), 'k', freq, db(Gbla), 'r', freq, db(varG1G2)/2, 'k:', freq, db(varGbla)/2, 'r:')
axis([0 0.5 -20 80])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(132)
plot(freq, (angle(Gbla1.*Gbla2))*180/pi, 'k', freq, (angle(Gbla))*180/pi, 'r')
axis([0 0.5 -200 200])
xlabel('f/fs')
ylabel('arg(G_{BLA}) (o)')

subplot(133)
plot(freq, db(Gbla-Gbla1.*Gbla2),'k', freq, db(varGbla + varG1G2)/2, 'r')
axis([0 0.5 -20 80])
xlabel('f/fs')
ylabel('error (dB)')

annotation('textbox',[0.35 0.8 0.05 0.2],'Color','k','String','Black: cascade BLA''s; Red: BLA cascade',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

zoom on;
shg

figure(2)
set(gcf, 'Position',[50 500 1200 300])              % larger window for plotting
subplot(131)
coh2 = 1./(1 + M*varGbla./abs(Gbla.^2));
plot(freq, db(Gbla1.*Gbla2), 'k', freq, db(Gbla), 'r', freq, db(coh2)/2, 'k--')
axis([0 0.5 -40 80])
xlabel('f/fs')
ylabel('amplitude (dB)')
title('Black: cascade BLA''s; red: BLA cascade; black dashes: coherence function')

subplot(132)
plot(freq, db(coh2)/2, 'k--')
axis([0.17 0.23 -10 0])
xlabel('f/fs')
ylabel('coherence (dB)')

subplot(133)
plot(freq, db(Gbla-Gbla1.*Gbla2),'k', freq, db(varGbla + varG1G2)/2, 'r')
axis([0.17 0.23 0 30])
xlabel('f/fs')
ylabel('error (dB)')
