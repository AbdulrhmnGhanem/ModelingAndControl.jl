% Chapter 5 Exercise 87.b and 87.c
% Exercise 87.b: Properties of output residuals - dynamic NL system
%          87.c: Predictive power of BLA - dynamic NL system
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                    %% 
% Predictive power of the best linear approximation  %%
% of a generalised Wiener-Hammerstein system, where  %%
% u is an odd random phase multisine                 %%
%                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                 % sampling frequency
N = 512*2;                                              % number of data points in one period (multisine) or block (Gaussian noise)
M = 10000;                                              % number of realisations (multisine) or blocks (Gaussian noise)
Mvalid = 1000;                                          % number of realisations for the validation

% definition first LTI system
a1 = [1 -0.5 0.9];
b1 = 1;

% definition second LTI system
a2 = [1 -1.5 0.7];
b2 = [0 1 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd multisine
FreqIndex = [1:2:N/2-1].';                                              % multisine exciting all odd harmonics Nyquist not included
freq = FreqIndex/N;                                                     % excited frequencies
F = length(freq);                                                       % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq);                                           % z^-1

freq_all = ([0:1:N-1]/N).';                                             % frequencies over complete unit circle
q_all = exp(-sqrt(-1)*2*pi*freq_all);                                   % z^-1 over the complete unit circle

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system
G2_all = polyval(fliplr(b2), q_all)./polyval(fliplr(a2), q_all);        % FRF second LTI system over complete unit circle

Gall = zeros(M, F);
Yall = zeros(M, N);
Uall = zeros(M, N);
Ysall = zeros(M, N);

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = exp(sqrt(-1)*2*pi*rand(size(freq)));               % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);
    U = fft(u)/sqrt(N);
    Uall(ii, :) = U.';
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% static nonlinear system with even NL contribution
 	z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
    Z = fft(z)/sqrt(N);
        
	% second LTI system
    Y = zeros(N,1);
    Y = Z.*G2_all;                                                              % steady state response at the excited frequencies
    Yall(ii, :) = Y.';

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                           % FRF of ii th realisation at the excited frequencies
end

stdu = sqrt(mean(u.^2));                                                        % rms value u(t)
Gbla = (mean(Gall, 1)).';                                                       % estimate best linear approximation (BLA)
clear Gall


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation stochastic nonlinear distortions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ysall(:, 1:N/2+1) = Yall(:, 1:N/2+1);                                                                                       % at the non-excited lines Ys equals Y
Ysall(:, FreqIndex+1) = Yall(:, FreqIndex+1) - Uall(:, FreqIndex+1) .* repmat(Gbla.', [M, 1]);                              % at the excited lines Ys = Y - Gbla*U
varYS = mean(abs(Ysall(:, FreqIndex+1)).^2, 1);                                % variance YS(k) at the excited frequencies
ys_all = 2*real(ifft(Ysall, [], 2)) * sqrt(N);                                                                              % stochastic nonlinear distortions for all realisations
varys = std(ys_all, [], 1).^2;                                                  % variance ys(t) over the realisations
varys_std = std(varys);                                                         % standard deviation var(ys(t)) over time
varys_mean = mean(varys);                                                       % mean value var(ys(t)) over time

stdys = sqrt(varys_mean);                                                       % standard deviation
% fraction inside 95% interval
fracVarys = length(find(abs(varys-varys_mean) <= 2*varys_std))/length(varys)*100;

% pdf stochastic nonlinear distortions
[ee_pre, xx_pre] = hist(ys_all(:), sqrt(N*Mvalid)/10);                          % histogram of u(t)
pdf_res_bla = ee_pre/(sum(ee_pre)*(mean(diff(xx_pre))));                        % experimental pdf = histogram divided by its area
pdf_res_gauss = 1/sqrt(2*pi)/stdys*exp(-xx_pre.^2/(2*stdys^2));                 % gaussian pdf with same mean and variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample cross-correlation between the input u(t) %
% and residual ys(t):                             %
%                                                 %
%      1/N sum u(t) * ys(t-tau)                   %
%       t = 0:N-1                                 %
%                                                 %
% averaged over MM realisations.                  %
% Since u(t) and ys(t) are periodic signals the   %
% periodic cross-correlation is calculated via    %
% DTF's                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MM = 1000;                                                                       % number of realisations taken to calculate Ruys
RR = zeros(MM, N);
for ii = 1:MM
   RR(ii, :) =  ifft(Uall(ii,:) .* conj(fft(ys_all(ii,:))))/sqrt(N);
end % ii
Ruys = mean(RR, 1);                                                              % mean value Ruys over the MM realisations
stdRuys = std(RR, [], 1)/sqrt(MM);                                               % standard deviation of the mean value
fracRuys = length(find((abs(Ruys) <= 2*stdRuys)))/length(Ruys)*100;              % fraction inside the 95% confidence interval


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auto-correlation stochastic nonlinear     %
% distortions at the F excited frequencies: %
%                                           %
%    1/(F-k) sum Ys(l) * conj(Ys(l-k))      %
%           l = k:F                         %
%                                           %
% averaged over MM realisations.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ysall = Ysall(:, FreqIndex+1);

% mean auto-correlation last MM realisations
MaxLag = 100;
MM = 1000;                                                                       
RR = zeros(MM, 2*MaxLag+1);
for  ii = 1:MM
    [RYs, Lags] = xcorr(Ysall(ii, :), MaxLag, 'unbiased');
    RR(ii, :) = RYs;
end % ii
RYs = mean(RR, 1);                                                              % mean value auto-correlation RYs
stdRYs = std(RR, [], 1)/sqrt(MM);                                               % standard deviation mean value
p = 0.95;                                                                       % 95% confidence level for circular complex normally distributed noise
ConfRYs_95 = sqrt(-log(1-p)) * stdRYs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probability density function Ys(k) %
% at one particular frequency        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The_k = floor(length(FreqIndex)/2);
[his_ampl, kk_ampl] = hist(abs(Ysall(:, The_k)), sqrt(M));                       % histogram of |Ys(k)|
pdf_ampl = his_ampl/(sum(his_ampl)*(mean(diff(kk_ampl))));                       % experimental pdf = histogram divided by its area
[his_phase, kk_phase] = hist(angle(Ysall(:, The_k)), sqrt(M));                   % histogram of angle(Ys(k))
pdf_phase = his_phase/(sum(his_phase)*(mean(diff(kk_phase))));                   % experimental pdf = histogram divided by its area
stdYsall = std(Ysall, [], 1);                                                    % standard deviation stochastic NL at the excited frequencies
stdYs_his = stdYsall(The_k);
pdf_ampl_Rayl = 2/stdYs_his^2*exp(-kk_ampl.^2/stdYs_his^2) .* kk_ampl;           % Rayleigh pdf with same second order moment
pdf_phase_unif = ones(size(kk_phase))/(2*pi);                                    % uniform pdf in [-pi, pi)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prediction output for the following inputs %
% multisines:                                %
% 1. random phase multisine                  %                  
% 2. Schroeder multisine                     %
% 3. zero phase multisine (= pulse signal)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% random phase multisine
U = zeros(N,1);
U(FreqIndex+1) = exp(sqrt(-1)*2*pi*rand(size(FreqIndex)));                      % random phase multisine with prescribed amplitude spectrum
u = 2*real(ifft(U))*sqrt(N);
U = fft(u)/sqrt(N);
% first LTI system
X = zeros(N,1);
X(FreqIndex+1) = U(FreqIndex+1).*G1;                                            % steady state response at the excited frequencies
x = 2*real(ifft(X))*sqrt(N);
% static nonlinear system with even NL contribution
z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
Z = fft(z)/sqrt(N);
% second LTI system
Y = zeros(N,1);
Y = Z.*G2_all;                                                                  % steady state response at the excited frequencies
y = real(ifft(Y))*sqrt(N);
Ypred = zeros(size(Y));
Ypred(FreqIndex+1) = U(FreqIndex+1) .* Gbla;
ypred = 2*real(ifft(Ypred))*sqrt(N);
erms = sqrt(mean((y-ypred).^2));                                                % rms value prediction error

% Schroeder phase multisine
U = zeros(N,1);
U(FreqIndex+1) = exp(sqrt(-1)*pi*FreqIndex.*(FreqIndex-1)/(FreqIndex(end)));   % Schroeder phase multisine with prescribed amplitude spectrum
u = 2*real(ifft(U))*sqrt(N);
U = fft(u)/sqrt(N);
% first LTI system
X = zeros(N,1);
X(FreqIndex+1) = U(FreqIndex+1).*G1;                                            % steady state response at the excited frequencies
x = 2*real(ifft(X))*sqrt(N);
% static nonlinear system with even NL contribution
z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
Z = fft(z)/sqrt(N);
% second LTI system
Y = zeros(N,1);
Y = Z.*G2_all;                                                                  % steady state response at the excited frequencies
y_Schr = real(ifft(Y))*sqrt(N);
Ypred = zeros(size(Y));
Ypred(FreqIndex+1) = U(FreqIndex+1) .* Gbla;
ypred_Schr = 2*real(ifft(Ypred))*sqrt(N);
erms_Schr = sqrt(mean((y_Schr-ypred_Schr).^2));                                 % rms value prediction error

% zero phase multisine (= pulse)
U = zeros(N,1);
U(FreqIndex+1) = ones(size(FreqIndex));                                         % zero phase multisine with prescribed amplitude spectrum
u = 2*real(ifft(U))*sqrt(N);
U = fft(u)/sqrt(N);
% first LTI system
X = zeros(N,1);
X(FreqIndex+1) = U(FreqIndex+1).*G1;                                            % steady state response at the excited frequencies
x = 2*real(ifft(X))*sqrt(N);
% static nonlinear system with even NL contribution
z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
Z = fft(z)/sqrt(N);
% second LTI system
Y = zeros(N,1);
Y = Z.*G2_all;                                                                  % steady state response at the excited frequencies
y_pulse = real(ifft(Y))*sqrt(N);
Ypred = zeros(size(Y));
Ypred(FreqIndex+1) = U(FreqIndex+1) .* Gbla;
ypred_pulse = 2*real(ifft(Ypred))*sqrt(N);
erms_pulse = sqrt(mean((y_pulse-ypred_pulse).^2));                              % rms value prediction error


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variance YS(k) for smaller number of frequencies: N/2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN = N/2;

% odd multisine
FreqIndex = [1:2:NN/2-1].';                                             % multisine exciting all odd harmonics Nyquist not included
freq_N2 = FreqIndex/NN;                                                 % excited frequencies
F = length(freq_N2);                                                    % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_N2);                                        % z^-1

freq_all = ([0:1:NN-1]/NN).';                                           % frequencies over complete unit circle
q_all = exp(-sqrt(-1)*2*pi*freq_all);                                   % z^-1 over the complete unit circle

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system
G2_all = polyval(fliplr(b2), q_all)./polyval(fliplr(a2), q_all);        % FRF second LTI system over complete unit circle

Gall = zeros(M, F);
Yall = zeros(M, NN);
Uall = zeros(M, NN);
Ysall = zeros(M, NN);

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(NN,1);
    U(FreqIndex+1) = exp(sqrt(-1)*2*pi*rand(size(freq_N2)));            % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(NN);
    U = fft(u)/sqrt(NN);
    Uall(ii, :) = U.';
    
    % first LTI system
    X = zeros(NN,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(NN);
    
	% static nonlinear system with even NL contribution
 	z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
    Z = fft(z)/sqrt(NN);
        
	% second LTI system
    Y = zeros(NN,1);
    Y = Z.*G2_all;                                                              % steady state response at the excited frequencies
    Yall(ii, :) = Y.';

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                           % FRF of ii th realisation at the excited frequencies
end

stdu = sqrt(mean(u.^2));                                                        % rms value u(t)
Gbla = (mean(Gall, 1)).';                                                       % estimate best linear approximation (BLA)
clear Gall

% calculation stochastic nonlinear distortions  %

Ysall(:, 1:NN/2+1) = Yall(:, 1:NN/2+1);                                                                                       % at the non-excited lines Ys equals Y
Ysall(:, FreqIndex+1) = Yall(:, FreqIndex+1) - Uall(:, FreqIndex+1) .* repmat(Gbla.', [M, 1]);                              % at the excited lines Ys = Y - Gbla*U
varYS_N2 = mean(abs(Ysall(:, FreqIndex+1)).^2, 1);                                % variance YS(k) at the excited frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variance YS(k) for smaller number of frequencies: N/4 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN = N/4;

% odd multisine
FreqIndex = [1:2:NN/2-1].';                                              % multisine exciting all odd harmonics Nyquist not included
freq_N4 = FreqIndex/NN;                                                  % excited frequencies
F = length(freq_N4);                                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_N4);                                         % z^-1

freq_all = ([0:1:NN-1]/NN).';                                            % frequencies over complete unit circle
q_all = exp(-sqrt(-1)*2*pi*freq_all);                                    % z^-1 over the complete unit circle

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);                    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);                    % FRF second LTI system
G2_all = polyval(fliplr(b2), q_all)./polyval(fliplr(a2), q_all);        % FRF second LTI system over complete unit circle

Gall = zeros(M, F);
Yall = zeros(M, NN);
Uall = zeros(M, NN);
Ysall = zeros(M, NN);

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(NN,1);
    U(FreqIndex+1) = exp(sqrt(-1)*2*pi*rand(size(freq_N4)));            % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(NN);
    U = fft(u)/sqrt(NN);
    Uall(ii, :) = U.';
    
    % first LTI system
    X = zeros(NN,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(NN);
    
	% static nonlinear system with even NL contribution
 	z = x + 0.01 * tanh(x) .* [x(end); x(1:end-1)].^2;                                           
    Z = fft(z)/sqrt(NN);
        
	% second LTI system
    Y = zeros(NN,1);
    Y = Z.*G2_all;                                                              % steady state response at the excited frequencies
    Yall(ii, :) = Y.';

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                           % FRF of ii th realisation at the excited frequencies
end

stdu = sqrt(mean(u.^2));                                                        % rms value u(t)
Gbla = (mean(Gall, 1)).';                                                       % estimate best linear approximation (BLA)
clear Gall

% calculation stochastic nonlinear distortions  %

Ysall(:, 1:NN/2+1) = Yall(:, 1:NN/2+1);                                                                                       % at the non-excited lines Ys equals Y
Ysall(:, FreqIndex+1) = Yall(:, FreqIndex+1) - Uall(:, FreqIndex+1) .* repmat(Gbla.', [M, 1]);                              % at the excited lines Ys = Y - Gbla*U
varYS_N4 = mean(abs(Ysall(:, FreqIndex+1)).^2, 1);                                % variance YS(k) at the excited frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the fractions and rms values prediction error %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fracVarys
fracRuys

erms
erms_Schr
erms_pulse


%%%%%%%%%%%%%%%%%%%%
% Plot the figures %
%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[150 300 1000 500])              % larger window for plotting

subplot(331);
Samples = [0:N-1].';
plot(Samples, varys, 'k', Samples, ones(size(varys.'))*varys_mean, 'g', Samples, ones(size(varys.'))*(varys_mean+2*varys_std), 'r', ...
     Samples, ones(size(varys.'))*(varys_mean-2*varys_std), 'r');
axis([-50, 1070, 5.7e-3, 7.8e-3])
xlabel('t (samples)')
ylabel('var(y_s(t))')

subplot(332)
plot(xx_pre, pdf_res_bla, 'k', xx_pre, pdf_res_gauss, 'r')
axis([-0.5, 0.5, 0, 10])
xlabel('y_s')
ylabel('pdf(y_s)')

subplot(333)
plot(Samples, 1e4*Ruys, 'k', Samples, 1e4*2*stdRuys, 'r', Samples, -1e4*2*stdRuys, 'r')
axis([-50, 1070, -2.4, 2.4])
xlabel('\tau (samples)')
ylabel('10^4 x R_{uy_s}(\tau)')

subplot(334)
semilogy(Lags, abs(RYs), 'k', Lags, ConfRYs_95, 'r');
axis([-100, 100, 0, 1])
xlabel('k (lag)')
ylabel('|R_{Y_s}(k)|')
 
subplot(335)
plot(kk_ampl, pdf_ampl, 'k', kk_ampl, pdf_ampl_Rayl, 'r');
axis([0, 0.08, 0, 50])
xlabel('|Y_s|')
ylabel('pdf(|Y_s|)')
 
subplot(336)
plot(kk_phase, pdf_phase, 'k', kk_phase, pdf_phase_unif, 'r');
axis([-4, 4, 0, 0.2])
xlabel('\angleY_s')
ylabel('pdf(\angleY_s)')
 
subplot(337)
plot(freq_N4, db(varYS_N4)/2, 'k')
axis([0, 0.5, -50, 0])
xlabel('f/f_s')
ylabel('var(Y_S) (dB)')
 
subplot(338)
plot(freq_N2, db(varYS_N2)/2, 'k')
axis([0, 0.5, -50, 0])
xlabel('f/f_s')
ylabel('var(Y_S) (dB)')
 
subplot(339)
plot(freq, db(varYS)/2, 'k')
axis([0, 0.5, -50, 0])
xlabel('f/f_s')
ylabel('var(Y_S) (dB)')

zoom on;
shg

figure(2)
set(gcf, 'Position',[150 300 1000 500])              % larger window for plotting

% prediction random phase multisine
subplot(231)
plot(Samples, y, 'k', Samples, ypred, 'g', Samples, y-ypred,'r');
axis([-50, 1070, -12, 12])
xlabel('t (samples)')
ylabel('y(t)')
title('Random phase')

% prediction Schroeder phase multisine
subplot(232)
plot(Samples, y_Schr, 'k', Samples, ypred_Schr, 'g', Samples, y_Schr-ypred_Schr,'r');
axis([-50, 1070, -12, 12])
xlabel('t (samples)')
ylabel('y(t)')
title('Schroeder phase')

% prediction zero phase multisine (= pulse)
subplot(233)
plot(Samples, y_pulse, 'k', Samples, ypred_pulse, 'g', Samples, y_pulse-ypred_pulse,'r');
axis([-50, 1070, -50, 50])
xlabel('t (samples)')
ylabel('y(t)')
title('Zero phase')

% prediction random phase multisine
subplot(234)
plot(Samples, y-ypred,'r');
axis([-50, 1070, -1, 1])
xlabel('t (samples)')
ylabel('y(t)')

% prediction Schroeder phase multisine
subplot(235)
plot(Samples, y_Schr-ypred_Schr,'r');
axis([-50, 1070, -1, 1])
xlabel('t (samples)')
ylabel('y(t)')

% prediction zero phase multisine (= pulse)
subplot(236)
plot(Samples, y_pulse-ypred_pulse,'r');
axis([-50, 1070, -5, 5])
xlabel('t (samples)')
ylabel('y(t)')

zoom on
shg
