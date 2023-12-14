% Chapter 5 Exercise 83.a
% Influence of rms value and pdf on the BLA
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influence of the rms value and pdf on the best      %%
% linear approximation of a static nonlinear system   %%
%                                                     %%
% Illustrated on                                      %%
%   1. atan function                                  %%
%   2. function with dead zone                        %%
%                                                     %%
% for                                                 %%
%   1. zero mean white Gaussian noise                 %%
%   2. zero mean white uniform noise                  %%
%   3. zero mean sine distributed noise               %%
%                                                     %%
% Rik Pintelon and Johan Schoukens                    %% 
% March 22, 2007                                      %%
%                                                     %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%
% Input signals %
%%%%%%%%%%%%%%%%%

N = 400*10;                                             % number of data points
u = linspace(-3, 3, N).';                               % the points are linearly spaced in [-3, 3]
urms = [0.1, 1.8, 3];                                   % the different rms values
nrms = length(urms);                                    % number of rms values
ug = randn(N, 1) * urms;                                % zero mean white Gaussian noise
uu = ((rand(N, 1) - 0.5)*sqrt(12)) * urms;              % zero mean uniformly distributed noise
phi = rand(N, 1)*2*pi - pi;                             % argument sine wave uniformly distributed in [-pi, pi]
us = (sin(phi)*sqrt(2)) * urms;                         % sinewave with rms value urms


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA atan function for different rms values and different pdf's %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yatan = atan(u);                                                % the exact value
yatan_g = atan(ug);                                             % atan function for the gaussian inputs
yatan_u = atan(uu);                                             % atan function for the uniform inputs
yatan_s = atan(us);                                             % atan function for the sine inputs
BLA_atan_g = zeros(N, nrms);                                    % BLA for gaussian inputs calculated over the interval [-3, 3]
BLA_atan_u = zeros(N, nrms);                                    % BLA for uniform inputs calculated over the interval [-3, 3]
BLA_atan_s = zeros(N, nrms);                                    % BLA for sine inputs calculated over the interval [-3, 3]
Coeffatan_g = zeros(nrms, 1);
Coeffatan_u = zeros(nrms, 1);
Coeffatan_s = zeros(nrms, 1);
for ii = 1:nrms
    % gaussian noise
    Coeffatan_g(ii) = sum(yatan_g(:, ii).*ug(:,ii))./sum(ug(:,ii).^2);      % coefficients best linear approximation (BLA)
    BLA_atan_g(:, ii) = Coeffatan_g(ii)*u;                                  % value of the BLA over the interval [-3, 3]
    % uniform noise
    Coeffatan_u(ii) = sum(yatan_u(:, ii).*uu(:,ii))./sum(uu(:,ii).^2);      % coefficients best linear approximation (BLA)
    BLA_atan_u(:, ii) = Coeffatan_u(ii)*u;                                  % value of the BLA over the interval [-3, 3]
    % sinewave
    Coeffatan_s(ii,:) = sum(yatan_s(:, ii).*us(:,ii))./sum(us(:,ii).^2);    % coefficients best linear approximation (BLA)
    BLA_atan_s(:, ii) = Coeffatan_s(ii)*u;                                  % value of the BLA over the interval [-3, 3]
end % ii

% Plot the results
figure(1)
set(gcf, 'Position',[50 350 1200 400])
subplot(121)
plot(u, yatan,'k', u, BLA_atan_g, 'r', u, BLA_atan_u, 'b', u, BLA_atan_s, 'g')
axis([-4, 4, -2, 2])
xlabel('u');
ylabel('y');
title('atan', 'FontName','Helvetica LT Std', 'FontSize', 10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA dead zone function for different rms values and different pdf's %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb = 1;
% the exact value
Zero_u = find(abs(u) <= bb);
Pos_u = find(u > bb);
Neg_u = find(u < -bb);
ydz = u;
ydz(Zero_u) = 0;
ydz(Pos_u) = u(Pos_u) - bb;
ydz(Neg_u) = u(Neg_u) + bb;

% the exact value for Gaussian inputs
[Zero_row_ug, Zero_column_ug] = find(abs(ug) <= bb);
[Pos_row_ug, Pos_column_ug] = find(ug > bb);
[Neg_row_ug, Neg_column_ug] = find(ug < -bb);
ydz_g = ug;
ydz_g(Zero_row_ug, Zero_column_ug) = 0;
ydz_g(Pos_row_ug, Pos_column_ug) = ug(Pos_row_ug, Pos_column_ug) - bb;
ydz_g(Neg_row_ug, Neg_column_ug) = ug(Neg_row_ug, Neg_column_ug) + bb;

% the exact value for uniform inputs
[Zero_row_uu, Zero_column_uu] = find(abs(uu) <= bb);
[Pos_row_uu, Pos_column_uu] = find(uu > bb);
[Neg_row_uu, Neg_column_uu] = find(uu < -bb);
ydz_u = uu;
ydz_u(Zero_row_uu, Zero_column_uu) = 0;
ydz_u(Pos_row_uu, Pos_column_uu) = uu(Pos_row_uu, Pos_column_uu) - bb;
ydz_u(Neg_row_uu, Neg_column_uu) = uu(Neg_row_uu, Neg_column_uu) + bb;

% the exact value for sine inputs
[Zero_row_us, Zero_column_us] = find(abs(us) <= bb);
[Pos_row_us, Pos_column_us] = find(us > bb);
[Neg_row_us, Neg_column_us] = find(us < -bb);
ydz_s = us;
ydz_s(Zero_row_us, Zero_column_us) = 0;
ydz_s(Pos_row_us, Pos_column_us) = us(Pos_row_us, Pos_column_us) - bb;
ydz_s(Neg_row_us, Neg_column_us) = us(Neg_row_us, Neg_column_us) + bb;

BLA_dz_g = zeros(N, nrms);                                    % BLA for gaussian inputs calculated over the interval [-3, 3]
BLA_dz_u = zeros(N, nrms);                                    % BLA for uniform inputs calculated over the interval [-3, 3]
BLA_dz_s = zeros(N, nrms);                                    % BLA for sine inputs calculated over the interval [-3, 3]
Coeffdz_g = zeros(nrms, 1);
Coeffdz_u = zeros(nrms, 1);
Coeffdz_s = zeros(nrms, 1);
for ii = 1:nrms
    % gaussian noise
    Coeffdz_g(ii,:) = sum(ydz_g(:, ii).*ug(:,ii))./sum(ug(:,ii).^2);        % coefficients best linear approximation (BLA)
    BLA_dz_g(:, ii) = Coeffdz_g(ii)*u;                                      % value of the BLA over the interval [-3, 3]
    % uniform noise
    Coeffdz_u(ii,:) = sum(ydz_u(:, ii).*uu(:,ii))./sum(uu(:,ii).^2);        % coefficients best linear approximation (BLA)
    BLA_dz_u(:, ii) = Coeffdz_u(ii)*u;                                      % value of the BLA over the interval [-3, 3]
    % sinewave
    Coeffdz_s(ii,:) = sum(ydz_s(:, ii).*us(:,ii))./sum(us(:,ii).^2);        % coefficients best linear approximation (BLA)
    BLA_dz_s(:, ii) = Coeffdz_s(ii)*u;                                      % value of the BLA over the interval [-3, 3]
end % ii

% Plot the results
subplot(122)
plot(u, ydz,'k', u, BLA_dz_g, 'r', u, BLA_dz_u, 'b', u, BLA_dz_s, 'g')
axis([-4, 4, -2, 2])
xlabel('u');
ylabel('y');
title('dead zone', 'FontName','Helvetica LT Std', 'FontSize', 10);
