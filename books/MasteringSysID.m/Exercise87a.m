% Chapter 5 Exercise 87a
% Predictive power BLA - static NL system
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
% of the static nonlinear system y = u^3, where      %%
% u is zero mean white uniformly distributed noise   %%
%                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zero mean uniformly distributed input %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10000;
Npred = 10000;
% generate N+Npred points; first N for the estimation of the BLA;
% next Npred for prediction output using BLA
e = (rand(1, N+Npred+1) - 0.5)*2;               % zero mean, standard deviation = 1/sqrt(3)
u = e;                                          % zero mean white noise with standard deviation stduf
stdu = 1/sqrt(3);

%%%%%%%
% BLA %
%%%%%%%

% output static nonlinear function
y = u.^3;

% best linear approximation (y and u have zero mean)
% calculated on the first N points
SelectFirstN = [1:N];
Gbla = mean(y(SelectFirstN).*u(SelectFirstN))./mean(u(SelectFirstN).^2);
Gsmall = 0.9*Gbla;                              % gain smaller than BLA
Glarge = 1.1*Gbla;                              % gain larger than BLA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output residuals on complete record %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ys = y - Gbla * u;
stdys = std(ys(SelectFirstN), 1);                  % needed for cross-correlation test on first N samples

% smaller gain
ys_small = y - Gsmall * u;
stdys_small = std(ys_small(SelectFirstN), 1);      % needed for cross-correlation test on first N samples

% larger gain
ys_large = y - Glarge * u;
stdys_large = std(ys_large(SelectFirstN), 1);      % needed for cross-correlation test on first N samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-correlation output residuals with input %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MaxLag = 500;
[xcorr_ys_u, lags] = xcorr(ys(SelectFirstN), u(SelectFirstN), MaxLag, 'unbiased');
stdxcorr_bla = (stdys*stdu)./(length(u(SelectFirstN)) - lags).^0.5;
select = [MaxLag+2:2*MaxLag+1];
frac = length(find(abs(xcorr_ys_u(select))./(2*stdxcorr_bla(select)) <= 1))/length(select)*100          % check the fraction inside the 95% confidence bound

[xcorr_ys_small_u, lags] = xcorr(ys_small(SelectFirstN), u(SelectFirstN), MaxLag, 'unbiased');
stdxcorr_small = (stdys_small*stdu)./(length(u(SelectFirstN)) - lags).^0.5;

[xcorr_ys_large_u, lags] = xcorr(ys_large(SelectFirstN), u(SelectFirstN), MaxLag, 'unbiased');
stdxcorr_large = (stdys_large*stdu)./(length(u(SelectFirstN)) - lags).^0.5;


%%%%%%%%%%%%%%%%%%%%%
% Prediction output %
%%%%%%%%%%%%%%%%%%%%%

SelectNpred = [N+1:N+Npred];
ypred = Gbla * u(SelectNpred);
epred = y(SelectNpred) - ypred;
stdepred = std(epred, 1);
[stdys, stdepred]
% pdf residual
[uu, xx_bla] = hist(epred, sqrt(length(epred)));               % histogram of u(t)
pdf_res_bla = uu/(sum(uu)*(mean(diff(xx_bla))));                 % experimental pdf = histogram divided by its area


% small gain
ypred_small = Gsmall * u(SelectNpred);
epred_small = y(SelectNpred) - ypred_small;
stdepred_small = std(epred_small, 1);
[stdys_small, stdepred_small]
% pdf residual
[uu, xx_small] = hist(epred_small, sqrt(length(epred_small)));               % histogram of u(t)
pdf_res_small = uu/(sum(uu)*(mean(diff(xx_small))));                 % experimental pdf = histogram divided by its area

% large gain
ypred_large = Glarge * u(SelectNpred);
epred_large = y(SelectNpred) - ypred_large;
stdepred_large = std(epred_large, 1);
[stdys_large, stdepred_large]
% pdf residual
[uu, xx_large] = hist(epred_large, sqrt(length(epred_large)));               % histogram of u(t)
pdf_res_large = uu/(sum(uu)*(mean(diff(xx_large))));                 % experimental pdf = histogram divided by its area


%%%%%%%%%%%%
% Graphics %
%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[150 200 900 800])              % larger window for plotting

% cross-correlation input u(t) and output residual ys(t)
subplot(331)
plot(lags, 1000*xcorr_ys_u, 'k+-', lags, 2000*stdxcorr_bla, 'r',  lags, -2000*stdxcorr_bla, 'r')
axis([-1.1*MaxLag, 1.1*MaxLag, -3, 3]) 
xlabel('\tau (samples)')
ylabel('10^3 x R_{uy_s}(\tau)')
title('LA = BLA')

subplot(332)
plot(lags, 1000*xcorr_ys_small_u, 'k+-', lags, 2000*stdxcorr_small, 'r',  lags, -2000*stdxcorr_small, 'r')
axis([-1.1*MaxLag, 1.1*MaxLag, -5, 25]) 
xlabel('\tau (samples)')
ylabel('10^3 x R_{uy_s}(\tau)')
title('LA < BLA')

subplot(333)
plot(lags, 1000*xcorr_ys_large_u, 'k+-', lags, 2000*stdxcorr_large, 'r',  lags, -2000*stdxcorr_large, 'r')
axis([-1.1*MaxLag, 1.1*MaxLag, -25, 5]) 
xlabel('\tau (samples)')
ylabel('10^3 x R_{uy_s}(\tau)')
title('LA > BLA')

% output prediction
subplot(334)
plot(SelectNpred/1000, y(SelectNpred), 'k', SelectNpred/1000, ypred, 'r', SelectNpred/1000, epred, 'g')
axis([SelectNpred(1)/1000 - 1, SelectNpred(end)/1000 + 1, -1.1, 1.1])
xlabel('t (ksamples)')
ylabel('y(t)')

subplot(335)
plot(SelectNpred/1000, y(SelectNpred), 'k', SelectNpred/1000, ypred_small, 'r', SelectNpred/1000, epred_small, 'g')
axis([SelectNpred(1)/1000 - 1, SelectNpred(end)/1000 + 1, -1.1, 1.1])
xlabel('t (ksamples)')
ylabel('y(t)')

subplot(336)
plot(SelectNpred/1000, y(SelectNpred), 'k', SelectNpred/1000, ypred_large, 'r', SelectNpred/1000, epred_large, 'g')
axis([SelectNpred(1)/1000 - 1, SelectNpred(end)/1000 + 1, -1.1, 1.1])
xlabel('t (ksamples)')
ylabel('y(t)')

% probability density function (pdf) prediction error ys(t) 
subplot(337)
plot(xx_bla, pdf_res_bla,'k')
axis([-0.45, 0.45, 0, 10])
xlabel('y_s')
ylabel('pdf(y_s)')

subplot(338)
plot(xx_small, pdf_res_small,'k')
axis([-0.45, 0.45, 0, 10])
xlabel('y_s')
ylabel('pdf(y_s)')

subplot(339)
plot(xx_large, pdf_res_large,'k')
axis([-0.45, 0.45, 0, 10])
xlabel('y_s')
ylabel('pdf(y_s)')
