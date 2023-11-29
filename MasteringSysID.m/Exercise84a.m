function excercise84a
% Chapter 5 Exercise 84a
% Comparison of Gaussian noise and random phase multisine
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %%
% Comparison BLA obtained by Gaussian noise and multisine %%
% with the same power spectrum (colouring and rms value)  %%
%                                                         %%
% Considered system: Generalised Wiener-Hammerstein       %%
%                                                         %%
%       1. linear time invariant system                   %%
%       2. nonlinear FIR system                           %%
%                                                         %%
%                                                         %%
% Rik Pintelon and Johan Schoukens                        %% 
% April 11, 2007                                          %%
%                                                         %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using random phase multisine excitations %
% 1. uniformly distributed phases              %
% 2. phase 0 or pi with equal probability      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. uniformly distributed phases %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FreqIndex = [1:1:N/2-1].';                              % multisine exciting all harmonics Nyquist not included
freq_m = FreqIndex/N;                                   % excited frequencies
F = length(freq_m);                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_m);                         % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                     % matrix of the FRF's of all realisations

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);    % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*2*pi*rand(size(freq_m)));     % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);                                        % and same rms value as the filtered Gaussian noise excitation
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t)) * x(t-1)^2
	z = atan(x) .* [x(end);x(1:end-1)].^2;                              % since x is periodic x(t-1) is obtained via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                   % FRF of ii th realisation at the excited frequencies
end

G_m_u = (mean(Gall, 1)).';
varG_m_u = (std(Gall, 0, 1).^2/M).';                                     % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. phases 0 or pi with prob. 1/2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FreqIndex = [1:1:N/2-1].';                              % multisine exciting all harmonics Nyquist not included
freq_m = FreqIndex/N;                                   % excited frequencies
F = length(freq_m);                                     % number of excited frequencies 
q = exp(-sqrt(-1)*2*pi*freq_m);                         % z^-1
L = polyval(fliplr(cL), q)./polyval(fliplr(dL), q);     % filter transfer function for shaping the amplitude spectrum of the multisine
Gall = zeros(M, F);                                     % matrix of the FRF's of all realisations

% first LTI system
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);    % FRF first LTI system

% second LTI system
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);    % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = abs(L).*exp(sqrt(-1)*pi*round(rand(size(freq_m))));    % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);                                            % and same rms value as the filtered Gaussian noise excitation
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1;                                    % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% NFIR system z(t) = atan(x(t)) * x(t-1)^2
	z = atan(x) .* [x(end);x(1:end-1)].^2;                                  % since x is periodic x(t-1) is obtained via a circular shift
    Z = fft(z)/sqrt(N);
    
	% second LTI system
    Y = zeros(N,1);
    Y(FreqIndex+1) = Z(FreqIndex+1).*G2;                                    % steady state response at the excited frequencies

    Gall(ii, :) = (Y(FreqIndex+1)./U(FreqIndex+1)).';                       % FRF of ii th realisation at the excited frequencies
end

G_m_b = (mean(Gall, 1)).';
varG_m_b = (std(Gall, 0, 1).^2/M).';                                    % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA using Gaussian noise excitation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = randn(N*M+1, 1);                                                    % add one sample because one sample is lost for calcuating the
                                                                        % response of the NFIR system
u = filter(cL, dL, e);                                                  % Gaussian noise with same colouring power spectrum as multisine

% first LTI system
a1 = [1 -0.5 0.9];
x = filter(b1, a1, u);

% NFIR system z(t) = atan(x(t)) * x(t-1)^2
z = atan(x(2:end)) .* x(1:end-1).^2;                                    % one sample is lost for calculating x(t-1)

% second LTI system
y = filter(b2, a2, z);
u = u(2:end);                                                           % remove one sample to compensate for the shift of the output
                                                                        % over one sample w.r.t. the input

[G_g, varG_g, freq_g] = BLA_arb(y, u, N, fs);                               % FRF estimate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the different BLA's %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[200 300 450 600])              % larger window for plotting
subplot(311)
plot(freq_m, db(G_m_u), 'k', freq_m, db(varG_m_u)/2, 'k:', freq_g, db(G_g), 'r', freq_g, db(varG_g)/2, 'r:')
axis([0 0.5 -60 40])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(312)
plot(freq_m, unwrap(angle(G_m_u))*180/pi, 'k', freq_g, unwrap(angle(G_g))*180/pi, 'r')
axis([0 0.5 -350 0])
xlabel('f/fs')
xlabel('f/fs')
ylabel('arg(G_{BLA}) (o)')

subplot(313)
% interpolation multisine measurements to have it at the same
% frequency grid as the Gaussian noise measurements
G_m_u_interp = interp1(freq_m, G_m_u, freq_g);
varG_m_u_interp = interp1(freq_m, varG_m_u, freq_g);
plot(freq_g, db(varG_m_u_interp + varG_g)/2, 'k', freq_g, db(G_m_u_interp-G_g), 'r');       % variance(Gm-Gg) = var(Gm)+var(Gg)
axis([0 0.5 -60 20])
xlabel('f/fs')
ylabel('error (dB)')
title('Black: std. FRF; red: complex difference FRF''s')

annotation('textbox',[0.25 0.8 0.05 0.2],'Color','k','String','Black: multisine; Red: Gaussian noise',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

zoom on;
shg


figure(2)
set(gcf, 'Position',[700 300 450 600])              % larger window for plotting
subplot(311)
plot(freq_m, db(G_m_u), 'k', freq_m, db(varG_m_u)/2, 'k:', freq_m, db(G_m_b), 'r', freq_m, db(varG_m_b)/2, 'r:')
axis([0 0.5 -60 40])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(312)
plot(freq_m, unwrap(angle(G_m_u))*180/pi, 'k', freq_m, unwrap(angle(G_m_b))*180/pi, 'r')
axis([0 0.5 -350 0])
xlabel('f/fs')
ylabel('arg(G_{BLA}) (o)')

subplot(313)
plot(freq_m, db(varG_m_u + varG_m_b)/2, 'k', freq_m, db(G_m_u-G_m_b), 'r');                  % variance(Gm_u-Gm_b) = var(Gm_u)+var(Gm_b)
axis([0 0.5 -60 20])
xlabel('f/fs')
ylabel('error (dB)')
title('Black: std. FRF; red: complex difference FRF''s')

annotation('textbox',[0.1 0.8 0.05 0.2],'Color','k','String','Multisines; black: uniform phases; Red: 0/pi phases',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

zoom on;
shg


function [G, varG, freq, Cov_vecG] = BLA_arb(y, u, N, fs)
%
%   function [G, varG, freq, Cov_vecG] = BLA_arb(y, u, N, fs);
%
%       calculation best linear approximation via 
%       division cross-power spectrum by auto-power spectrum
%
%   Output
%       G           =   estimated frequency response function (FRF)
%       varG        =   variance estimated FRF
%       freq        =   frequencies at which FRF is estimated
%       Cov_vecG	=   covariance matrix of vec(G) (for MIMO systems only)
%
%   Input
%       y           =   output signal; size N*M x 1 or 1x N*M or ny x N*M
%       u           =   input signal; size N*M x 1 or 1 x N*M or nu x N*M
%       N           =   block length for FRF estimate
%       fs          =   sampling frequency
%
%   Rik Pintelon and Johan Schoukens
%   March 23, 2007
%   version 18 June 2010
%

% number of inputs and outputs
ny = min(size(y));
nu = min(size(u));

% number of blocks
M = floor(length(u)/N);

% frequencies at which the FRF is calculated
freq = ([1:1:N/2-2].' + 0.5)*fs/N;


if ny*nu == 1
    
    %%%%%%%%
    % SISO %
    %%%%%%%%
    
    u = u(:);
    y = y(:);

    u = reshape(u(1:M*N), N, M);
    y = reshape(y(1:M*N), N, M);

    u = fft(u, [], 1)/sqrt(N);
    u = u(2:N/2, :);
    y = fft(y, [], 1)/sqrt(N);
    y = y(2:N/2, :);

    % leakage suppression by taking the difference
    u = diff(u, 1, 1);
    y = diff(y, 1, 1);

    % FRF estimate Syu/Suu
    Syu = mean(y.*conj(u), 2);
    Suu = mean(abs(u.^2), 2);
    G = Syu ./ Suu;

    % variance estimate FRF
    Syy = mean(abs(y.^2), 2);
    varG = 1/(M-1) * (Syy - abs(Syu.^2)./Suu)./Suu;
    
else
    
    %%%%%%%%
    % MIMO %
    %%%%%%%%
    
    % put the vectors in the correct format
    if ny == 1
        y = y(:).';
    end % if ny = 1
    if nu == 1
        u = u(:).';
    end % if nu = 1
   
    % take fft of each segment
    y = y(:,1:M*N);
    u = u(:,1:M*N);
    y = reshape(y, [ny, N, M]);
    u = reshape(u, [nu, N, M]);
    y = permute(y, [1, 3, 2]);
    u = permute(u, [1, 3, 2]);
    Y = fft(y, [], 3)/sqrt(N);		% ny x M x N matrix
    U = fft(u, [], 3)/sqrt(N);
    Y = Y(:,:,2:N/2);
    U = U(:,:,2:N/2);

    % differentiate to suppress leakage
    Y = diff(Y, 1, 3);
    U = diff(U, 1, 3);
    F = size(Y, 3);

    % estimate output power spectrum = sum true output power + noise power
    Syy = zeros(ny, ny, F);
    for ii = 1:ny
        for jj = 1:ii
            Syy(ii, jj, :) = mean(Y(ii, :, :).*conj(Y(jj, :, :)), 2);
            Syy(jj, ii, :) = conj(Syy(ii, jj, :));
        end % jj
    end % ii

    % estimate input power spectrum
    Suu = zeros(nu, nu, F);
    for ii = 1:nu
        for jj = 1:ii
            Suu(ii, jj, :) = mean(U(ii, :, :).*conj(U(jj, :, :)), 2);
            Suu(jj, ii, :) = conj(Suu(ii, jj, :));
        end % jj
    end % ii

    % estimate cross power spectrum
    Syu = zeros(ny, nu, F);
    for ii = 1:ny
        for jj = 1:nu
            Syu(ii, jj, :) = mean(Y(ii, :, :).*conj(U(jj, :, :)), 2);
        end % jj
    end % ii

    % estimate output noise power spectrum and FRF
    Cy = zeros(ny, ny, F);
    G = zeros(ny, nu, F);
    Cov_vecG = zeros(ny*nu, ny*nu, N-1);
    for kk = 1:F
        invSuu = inv(squeeze(Suu(:, :, kk)));
        G(:, :, kk) = Syu(:, :, kk)*invSuu;
        Cy(:, :, kk) = Syy(:, :, kk) - G(:, :, kk)*Syu(:, :, kk)';
        % remove imaginary part on diagonal
        Cy(:, :, kk) = Cy(:, :, kk) - diag(sqrt(-1)*imag(diag(Cy(:, :, kk))));
        % cov(vec(G)) = kron(inv(Suu).', Cy) / M
        Cov_vecG(:, :, kk) = kron(invSuu.', squeeze(Cy(:, :, kk)));
    end % kk

    % correct for bias;not for diff because Suu and Cy are both 2 times
    % too large => factor two cancels in the ratio: 
    %       kron(invSuu.', squeeze(Cy(:, :, kk)))
    Cov_vecG = Cov_vecG/(M-nu);
    
    % variance entries FRF
    varG = zeros(ny, nu, F);
    for kk = 1:F
        varG(:, :, kk) = reshape(diag(squeeze(Cov_vecG(:, :, kk))), [ny, nu] );
    end % kk

end % if



