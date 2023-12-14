function Exercise83b
% Chapter 5 Exercise 83.b
% Influence of power spectrum coloring and pdf on the BLA
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %%
% Influence of the power spectrum and pdf on the best  %%
% linear approximation of y = u^3                      %%
%                                                      %%
% Considered inputs                                    %%
%                                                      %%
%       1. white noise                                 %%
%       2. filtered white noise                        %%
%                                                      %%
% Considered pdf's                                     %%
%                                                      %%
%       1. uniform                                     %%
%       2. normal                                      %%
%                                                      %%
% Rik Pintelon and Johan Schoukens                     %% 
% March 23, 2007                                       %%
%                                                      %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% general settings
fs = 1;                                     % sampling frequency
N = 512;                                    % number of points in one block
M = 8000;                                   % number of blocks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA for (filtered) uniform white noise %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = 2*rand(N*M,1) - 1;                      % uniform white noise [-1, 1]

% example 1.1 PhD thesis Martin Enqvist: minimum phase filtered uniform white noise
u = filter([1 0.5], [1], e);                % input signal
std(u)
y = u.^3;
[G1_u, varG1_u, freq] = BLA_arb(y, u, N, fs);   % FRF estimate

% example 1.2 PhD thesis Martin Enqvist: minimum phase filtered uniform white noise
u = filter([0.5 1], [1], e);                % input signal
std(u)
y = u.^3;
[G2_u, varG2_u, freq] = BLA_arb(y, u, N, fs);   % FRF estimate

% Gaussian white noise input
u = sqrt(5/4)*e;                            % uniform white noise with same std as the filtered white noise
std(u)
y = u.^3;
[G3_u, varG3_u, freq] = BLA_arb(y, u, N, fs);   % FRF estimate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA for (filtered) gaussian white noise %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = randn(N*M,1)/sqrt(3);                   % gaussian white noise with same std as in the uniform case

% example 1.1 PhD thesis Martin Enqvist: minimum phase filtered uniform white noise
u = filter([1 0.5], [1], e);                % input signal
std(u)
y = u.^3;
[G1_g, varG1_g, freq] = BLA_arb(y, u, N, fs);   % FRF estimate

% example 1.2 PhD thesis Martin Enqvist: minimum phase filtered uniform white noise
u = filter([0.5 1], [1], e);                % input signal
std(u)
y = u.^3;
[G2_g, varG2_g, freq] = BLA_arb(y, u, N, fs);   % FRF estimate

% uniform white noise input
u = sqrt(5/4)*e;                            % uniform white noise with same std as the filtered white noise
std(u)
y = u.^3;
[G3_g, varG3_g, freq] = BLA_arb(y, u, N, fs);   % FRF estimate


figure(1)
set(gcf, 'Position',[300 200 600 800])      % larger window for plotting
subplot(321)
plot(freq, db(G1_u), 'k', freq, db(G2_u), 'r', freq, db(G3_u), 'b')
xlabel('f/fs');
ylabel('|G_{BLA}| (dB)')
axis1 = [0 0.5 -6 4];
axis(axis1);
subplot(323)
plot(freq, angle(G1_u)*180/pi, 'k', freq, angle(G2_u)*180/pi, 'r', freq, angle(G3_u)*180/pi, 'b')
xlabel('f/fs');
ylabel('arg(G_{BLA}) (o)')
axis2 = [0 0.5 -20 20];
axis(axis2);
subplot(325)
plot(freq, db(varG1_u)/2, 'k', freq, db(varG2_u)/2, 'r', freq, db(varG3_u)/2, 'b')
xlabel('f/fs');
ylabel('var(G_{BLA}) (dB)')
axis3 = [0 0.5 -50 -30];
axis(axis3);
annotation('textbox',[0.09 0.8 0.05 0.2],'Color','k','String','(Filtered) uniform white noise',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

subplot(322)
plot(freq, db(G1_g), 'k', freq, db(G2_g), 'r', freq, db(G3_g), 'b')
xlabel('f/fs');
ylabel('|G_{BLA}| (dB)')
axis(axis1);
subplot(324)
plot(freq, angle(G1_g)*180/pi, 'k', freq, angle(G2_g)*180/pi, 'r', freq, angle(G3_g)*180/pi, 'b')
xlabel('f/fs');
ylabel('arg(G_{BLA}) (o)')
axis(axis2);
subplot(326)
plot(freq, db(varG1_g)/2, 'k', freq, db(varG2_g)/2, 'r', freq, db(varG3_g)/2, 'b')
xlabel('f/fs');
ylabel('var(G_{BLA}) (dB)')
axis(axis3);
annotation('textbox',[0.53 0.8 0.05 0.2],'Color','k','String','(Filtered) Gaussian white noise',...
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



