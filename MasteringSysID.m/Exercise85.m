function Exercise85
% Chapter 5 Exercise 85
% Influence of Even and Odd Nonlinearities on BLA
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %%
% BLA of even and odd nonlinear systems excited by odd   %%
% random phase multisines and chi2 distributed noise     %%
%                                                        %%
% Goal: influence even and odd nonlinearities on BLA and %%
%       its uncertainty                                  %%
%                                                        %%
% Rik Pintelon and Johan Schoukens                       %% 
% April 18, 2007                                         %%
% modified May 16, 2007                                  %%
%                                                        %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                                 % sampling frequency
N = 512*2;                                              % number of data points in one period (multisine) or block (Gaussian noise)
M = 10000;                                              % number of realisations (multisine) or blocks (Gaussian noise)

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
freq_o = FreqIndex/N;                                                   % excited frequencies
F = length(freq_o);                                                     % number of excited frequencies 
q_o = exp(-sqrt(-1)*2*pi*freq_o);                                       % z^-1
Gall_E = zeros(M, F);                                                   % matrix of the FRF's of all realisations
Gall_O = zeros(M, F);                                                   % matrix of the FRF's of all realisations

% first LTI system
G1_o = polyval(fliplr(b1), q_o)./polyval(fliplr(a1), q_o);              % FRF first LTI system

% second LTI system
G2_o = polyval(fliplr(b2), q_o)./polyval(fliplr(a2), q_o);              % FRF second LTI system

% calculation FRF over the M random phase realisations multisine 
for ii=1:M
    home, ii
    U = zeros(N,1);
    U(FreqIndex+1) = exp(sqrt(-1)*2*pi*rand(size(freq_o)));             % random phase multisine with prescribed amplitude spectrum
    u = 2*real(ifft(U))*sqrt(N);
    U = fft(u)/sqrt(N);
    
    % first LTI system
    X = zeros(N,1);
	X(FreqIndex+1) = U(FreqIndex+1).*G1_o;                              % steady state response at the excited frequencies
    x = 2*real(ifft(X))*sqrt(N);
    
	% static nonlinear system with even NL contribution
 	z_E = x + abs(atan(x) .* [x(end);x(1:end-1)].^2);                                           
    Z_E = fft(z_E)/sqrt(N);
    
	% static nonlinear system with odd NL contribution
 	z_O = x + atan(x) .* [x(end);x(1:end-1)].^2;                                           
    Z_O = fft(z_O)/sqrt(N);
    
	% second LTI system
    Y_E = zeros(N,1);
    Y_E(FreqIndex+1) = Z_E(FreqIndex+1).*G2_o;                          % steady state response at the excited frequencies

	% second LTI system
    Y_O = zeros(N,1);
    Y_O(FreqIndex+1) = Z_O(FreqIndex+1).*G2_o;                          % steady state response at the excited frequencies

    Gall_E(ii, :) = (Y_E(FreqIndex+1)./U(FreqIndex+1)).';               % FRF of ii th realisation at the excited frequencies
    Gall_O(ii, :) = (Y_O(FreqIndex+1)./U(FreqIndex+1)).';               % FRF of ii th realisation at the excited frequencies
end

G_E_o = (mean(Gall_E, 1)).';
varG_E_o = (std(Gall_E, 0, 1).^2/M).';                                  % variance mean value

G_O_o = (mean(Gall_O, 1)).';
varG_O_o = (std(Gall_O, 0, 1).^2/M).';                                  % variance mean value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA's using chi2 distributed noise excitation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = (randn(N*M+1, 1).^2 - 1)/2;                                         % zero mean chi2 distributed noise with standard deviation 1/sqrt(2)
                                                                        % add one sample because one sample is lost for calcuating the
                                                                        % response of the NFIR system
% first LTI system
x = filter(b1, a1, u);

% static nonlinear system with even NL contribution
z_E = x(2:end) + abs(atan(x(2:end)) .* x(1:end-1).^2);                  % one sample is lost for calculating x(t-1)                                   

% static nonlinear system with odd NL contribution
z_O = x(2:end) + atan(x(2:end)) .* x(1:end-1).^2;                       % one sample is lost for calculating x(t-1)                             

% second LTI system
y_E = filter(b2, a2, z_E);
y_O = filter(b2, a2, z_O);
u = u(2:end);                                                           % remove one sample to compensate for the shift of the output
                                                                        % over one sample w.r.t. the input

[G_E_chi2, varG_E_chi2, freq_chi2] = BLA_arb(y_E, u, N, fs);                % FRF estimate for even nonlinearity
[G_O_chi2, varG_O_chi2, freq_chi2] = BLA_arb(y_O, u, N, fs);                % FRF estimate for odd nonlinearity

q_chi2 = exp(-sqrt(-1)*2*pi*freq_chi2);                                 % z^-1
% first LTI system
G1_chi2 = polyval(fliplr(b1), q_chi2)./polyval(fliplr(a1), q_chi2);     % FRF first LTI system

% second LTI system
G2_chi2 = polyval(fliplr(b2), q_chi2)./polyval(fliplr(a2), q_chi2);     % FRF second LTI system


% pdf output first dynamic system
stdx = sqrt(mean(x.^2));                                                % mean value is zero => estimate standard deviation
[hx, xx] = hist(x(:), sqrt(length(x(:)))/10);                           % histogram of x(t)
pdf_x = hx/(sum(hx)*(mean(diff(xx))));                                  % experimental pdf = histogram divided by its area
pdf0_x = exp(-xx.^2/(2*stdx^2))/(sqrt(2*pi)*stdx);                      % gaussian pdf with zero mean and variance stdu^2

m1_u = mean(u)*sqrt(2);                                                 % mean relative to the standard deviation
m3_u = mean(u.^3)*sqrt(2)^3;                                            % third order moment relative to the standard deviation 
                                                                        % = measure of the skewness of u 
m1_x = mean(x)/stdx;                                                    % mean relative to the standard deviation
m3_x = mean(x.^3)/stdx^3;                                               % third order moment relative to the standard deviation 
                                                                        % = measure of the skewness of x
m3_g = mean(randn(size(u)).^3);                                         % sample third order moment normal distribution
disp(sprintf('Sample third order moment u(t) relative to the standard deviation: %4.3e', m3_u))
disp(sprintf('Sample third order moment x(t) relative to the standard deviation: %4.3e', m3_x))
disp(sprintf('Sample third order moment gaussian relative to the standard deviation: %4.3e', m3_g))


%%%%%%%%%%%%%%%%%%%%%
% Compare the BLA's %
%%%%%%%%%%%%%%%%%%%%%

% even nonlinear distortion; odd random phase multisine
figure(1)
set(gcf, 'Position',[50 300 1200 600])              % larger window for plotting
subplot(231)
plot(freq_o, db(G_E_o), 'k', freq_o, db(G1_o.*G2_o), 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(232)
plot(freq_o, unwrap(angle(G_E_o))*180/pi, 'k', freq_o, unwrap(angle(G1_o.*G2_o))*180/pi, 'r')
xlabel('f/fs')
ylabel('arg(G_{BLA}) (0)')

subplot(233)
plot(freq_o, db(G_E_o-G1_o.*G2_o), 'k')
axis([0 0.5 -300 -200])
xlabel('f/fs')
ylabel('error (dB)')

% odd nonlinear distortion; odd random phase multisine
subplot(234)
plot(freq_o, db(G_O_o), 'k', freq_o, db(G1_o.*G2_o), 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(235)
plot(freq_o, unwrap(angle(G_O_o))*180/pi, 'k', freq_o, unwrap(angle(G1_o.*G2_o))*180/pi, 'r')
axis([0 0.5 -300 0]);
xlabel('f/fs')
ylabel('arg(G_{BLA}) (0)')

subplot(236)
plot(freq_o, db(G_O_o-G1_o.*G2_o), 'k', freq_o, db(varG_O_o)/2, 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('error (dB)')

annotation('textbox',[0.15 0.8 0.05 0.2],'Color','k','String','Results odd multisine. Black: BLA; red: underlying linear system; First row: Even nonlinearity; Second row: odd nonlinearity',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

zoom on;
shg


% even nonlinear distortion; chi2 noise
figure(2)
set(gcf, 'Position',[50 300 1200 600])              % larger window for plotting
subplot(231)
plot(freq_chi2, db(G_E_chi2), 'k', freq_chi2, db(G1_chi2.*G2_chi2), 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(232)
plot(freq_chi2, unwrap(angle(G_E_chi2))*180/pi, 'k', freq_chi2, unwrap(angle(G1_chi2.*G2_chi2))*180/pi, 'r')
axis([0 0.5 -400 0]);
xlabel('f/fs')
ylabel('arg(G_{BLA}) (o)')

subplot(233)
plot(freq_chi2, db(G_E_chi2-G1_chi2.*G2_chi2), 'k', freq_chi2, db(varG_E_chi2)/2, 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('error (dB)')

% odd nonlinear distortion
subplot(234)
plot(freq_chi2, db(G_O_chi2), 'k', freq_chi2, db(G1_chi2.*G2_chi2), 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('G_{BLA} (dB)')

subplot(235)
plot(freq_chi2, unwrap(angle(G_O_chi2))*180/pi, 'k', freq_chi2, unwrap(angle(G1_chi2.*G2_chi2))*180/pi, 'r')
axis([0 0.5 -400 0]);
xlabel('f/fs')
ylabel('arg(G_{BLA}) (o)')

subplot(236)
plot(freq_chi2, db(G_O_chi2-G1_chi2.*G2_chi2), 'k', freq_chi2, db(varG_O_chi2)/2, 'r')
axis([0 0.5 -50 50])
xlabel('f/fs')
ylabel('error (dB)')

annotation('textbox',[0.15 0.8 0.05 0.2],'Color','k','String','Results Chi2 noise. Black: BLA; red: underlying linear system; First row: Even nonlinearity; Second row: odd nonlinearity',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

zoom on;
shg

% comparison pdf x with normal distribution with same mean value and
% variance
figure(3)
plot(xx, pdf_x, 'k', xx, pdf0_x,'r')
xlabel('x')
ylabel('pdf(x)')
title('Black: pdf(x); red: normal distribution')




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


