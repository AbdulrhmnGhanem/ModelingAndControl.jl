function exercise83c
% Chapter 5 Exercise 83c
% Influence of length of impulse response of signal filter on the BLA
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %%
% Influence of length of the impulse response of the   %%
% signal filter on the best linear approximation of    %%
%                                                      %%
%                      y = u^3                         %%
%                                                      %%
% The input is of the form                             %%
%                                                      %%
%                     u(t) = H(q) e(t)                 %%
%                                                      %%
% with e(t) white uniformly distributed noise          %%
%                                                      %%
%                                                      %%
% Rik Pintelon and Johan Schoukens                     %% 
% April 4, 2007                                        %%
%                                                      %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                             % sampling frequency
N = 512;                                            % number of points in one block
M = 8000;                                           % number of blocks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA for (filtered) uniform white noise %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = 2*rand(N*M,1) - 1;                              % uniform white noise [-1, 1]

% uniform white noise input
u = sqrt(5/4)*e;                                    % uniform white noise with same std as the filtered white noise
stdu = sqrt(5/4/3);                                 % standard deviation u(t)
[uu1, xx1] = hist(u, sqrt(length(u)));              % histogram of u(t)
pdf1 = uu1/(sum(uu1)*(mean(diff(xx1))));            % experimental pdf = histogram divided by its area
pdf10 = exp(-xx1.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);  % gaussian pdf with zero mean and variance stdu^2
y = u.^3;
[G1_u, varG1_u, freq] = BLA_arb(y, u, N, fs);           % FRF estimate

% filtered uniform white noise
fc = 0.25;                                          % cutoff frequency
[b2, a2] = cheby1(2, 1, 2*fc);                      % chebyshev filter of order 2, 1 dB passband ripple, cutoff frequency fc
u = filter(b2, a2, e);                              % input signal
u = u/std(u)*stdu;                                  % same rms value as the white noise case
[uu2, xx2] = hist(u, sqrt(length(u)));              % histogram of u(t)
pdf2 = uu2/(sum(uu2)*(mean(diff(xx2))));            % experimental pdf = histogram divided by its area
pdf20 = exp(-xx2.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);  % gaussian pdf with zero mean and variance stdu^2
y = u.^3;
[G2_u, varG2_u, freq] = BLA_arb(y, u, N, fs);           % FRF estimate
Select2 = find(freq <= fc);                         % select excited frequency band

% filtered uniform white noise
fc = 0.25;                                          % cutoff frequency
[b3, a3] = cheby1(35, 1, 2*fc);                     % chebyshev filter of order 35, 1 dB passband ripple, cutoff = fc
u = filter(b3, a3, e);                              % input signal
u = u/std(u)*stdu;                                  % same rms value as the white noise case
[uu3, xx3] = hist(u, sqrt(length(u)));              % histogram of u(t)
pdf3 = uu3/(sum(uu3)*(mean(diff(xx3))));            % experimental pdf = histogram divided by its area
pdf30 = exp(-xx3.^2/(2*stdu^2))/(sqrt(2*pi)*stdu);  % gaussian pdf with zero mean and variance stdu^2
y = u.^3;
[G3_u, varG3_u, freq] = BLA_arb(y, u, N, fs);           % FRF estimate
Select3 = find(freq <= fc);                         % select excited frequency band


figure(1)
set(gcf, 'Position',[300 100 300 800])              % larger window for plotting
subplot(311)
plot(freq, db(G1_u), 'k', freq(Select2), db(G2_u(Select2)), 'r', freq(Select3), db(G3_u(Select3)), 'b')
xlabel('f/fs');
ylabel('|G_{BLA}| (dB)')
axis1 = [0 0.5 -6 4];
axis(axis1);
subplot(312)
plot(freq, angle(G1_u)*180/pi, 'k', freq(Select2), angle(G2_u(Select2))*180/pi, 'r', freq(Select3), angle(G3_u(Select3))*180/pi, 'b')
xlabel('f/fs');
ylabel('arg(G_{BLA}) (o)')
axis2 = [0 0.5 -20 20];
axis(axis2);
subplot(313)
plot(freq, db(varG1_u)/2, 'k', freq(Select2), db(varG2_u(Select2))/2, 'r', freq(Select3), db(varG3_u(Select3))/2, 'b')
xlabel('f/fs');
ylabel('var(G_{BLA}) (dB)')
axis3 = [0 0.5 -50 -30];
axis(axis3);
annotation('textbox',[0.1 0.8 0.05 0.2],'Color','k','String','(Filtered) uniform white noise',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(2);
set(gcf, 'Position',[650 100 300 800])              % larger window for plotting
subplot(311)
plot(xx1(1:end-1), pdf1(1:end-1), xx1, pdf10);
axis([-3 3 0 0.8])
subplot(312)
plot(xx2, pdf2, xx2, pdf20);
axis([-3 3 0 0.8])
subplot(313)
plot(xx3, pdf3, xx3, pdf30);
axis([-3 3 0 0.8])
annotation('textbox',[0.1 0.8 0.05 0.2],'Color','k','String','pdf of filtered uniform white noise',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );
       
       
       
       
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



