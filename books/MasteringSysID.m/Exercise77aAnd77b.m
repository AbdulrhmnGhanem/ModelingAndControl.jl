% Chapter 5 Exercise 77.a and 77.b
% Exercise 77.a: Single sine response of a static nonlinear system
% Exercise 77.b: Multisine response of a static nonlinear system
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 6 December 2010



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %% 
% Response static nonlinear system to a periodic input %%
%                                                      %%
%   Systems considered                                 %%
%       1. u^2 and cosh(u)      (even nonlinearity)    %%
%       2. u^3 and sinh(u)      (odd nonlinearity)     %%
%       3. alfa*u^2 + beta*u^3                         %%
%                                                      %%
%   Inputs considered                                  %%
%       1. single sine                                 %%
%       2. full multisine                              %%
%       3. odd multisine                               %%
%                                                      %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% single sine response %
%%%%%%%%%%%%%%%%%%%%%%%%

A = 2;                                  % amplitude sinewave
f0 = 1;                                 % frequency sinewave
fs = 200;                               % sampling frequency
N = fs/f0;                              % number of points in one period
t = [0:N-1].';                          % time in samples
u = A*cos(2*pi*t/N);                    % single sine input

% even nonlinearities
y2 = u.^2;                              % square
ych = cosh(u);                          % hyperbolic cosine

% cubic nonlinearity
y3 = u.^3;                              % cubic
ysh = sinh(u);                          % hyperbolic sine

% calculation of the spectra
U = fft(u)/N;                           % scaling with N to recover the correct Fourier coefficient
Y2 = fft(y2)/N;
Ych = fft(ych)/N;
Y3 = fft(y3)/N;
Ysh = fft(ysh)/N;

% selection DFT lines shown in the figures
DFTlines = [0:20];
U = U(DFTlines+1);
Y2 = Y2(DFTlines+1);
Ych = Ych(DFTlines+1);
Y3 = Y3(DFTlines+1);
Ysh = Ysh(DFTlines+1);

% plot the results square
figure(1)
subplot(221)
plot(t, u);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y2);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(U), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y2), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Squaring Nonlinear Static System - Single Sine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );


% plot the results cosh
figure(2)
subplot(221)
plot(t, u);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, ych);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(U), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Ych), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Cosh Nonlinear Static System - Single Sine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );


% plot the results cubic
figure(3)
subplot(221)
plot(t, u);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y3);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(U), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y3), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Cubic Nonlinear Static System - Single Sine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );


% plot the results sinh
figure(4)
subplot(221)
plot(t, u);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, ysh);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(U), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Ysh), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Sinh Nonlinear Static System - Single Sine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

       
%%%%%%%%%%%%%%%%%%%%%%%
% multi sine response %
%%%%%%%%%%%%%%%%%%%%%%%

% full multisine with 5 components: [1:5] * f0
Uf = zeros(1, N).';
kkf = [1:5].';
Uf(kkf+1) = exp(-sqrt(-1)*pi*kkf.*(kkf-1)/length(kkf));
uf = real(ifft(Uf));
uf = uf/std(uf, 1);                            % rms value equal to one

% odd multisine with 3 components: [1, 5] * f0
Uo = zeros(1, N).';
kko = [1:2:5].';
Uo(kko+1) = exp(-sqrt(-1)*pi*kko.*(kko-1)/length(kko));
uo = real(ifft(Uo));
uo = uo/std(uo, 1);                            % rms value equal to one

% response static systems to full and odd multisines
y2f = uf.^2;
y2o = uo.^2;
y3f = uf.^3;
y3o = uo.^3;

% calculation of the spectra
Uf = fft(uf)/N;                               % scaling with N to get the Fourier coefficients
Uo = fft(uo)/N;
Y2f = fft(y2f)/N;
Y2o = fft(y2o)/N;
Y3f = fft(y3f)/N;
Y3o = fft(y3o)/N;

% selection DFT lines shown in the figures
eeps = 1e-15;                                 % if exactly zero then the spectral line is not shown
DFTlines = [0:20];
Uf = Uf(DFTlines+1) + eeps;
Uo = Uo(DFTlines+1) + eeps;
Y2f = Y2f(DFTlines+1) + eeps;
Y2o = Y2o(DFTlines+1) + eeps;
Y3f = Y3f(DFTlines+1) + eeps;
Y3o = Y3o(DFTlines+1) + eeps;

% plot the results square
figure(5)
subplot(221)
plot(t, uf);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y2f);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uf), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y2f), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Squaring Nonlinear Static System - Full Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(6)
subplot(221)
plot(t, uo);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y2o);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uo), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y2o), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Squaring Nonlinear Static System - Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );


% plot the results cubic
figure(7)
subplot(221)
plot(t, uf);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y3f);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uf), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y3f), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Cubic Nonlinear Static System - Full Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(8)
subplot(221)
plot(t, uo);
axis([-5, 205, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, y3o);
axis([-5, 205, -11, 11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uo), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(Y3o), 'k+')
axis([-0.5, 20.5, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Cubic Nonlinear Static System - Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );
