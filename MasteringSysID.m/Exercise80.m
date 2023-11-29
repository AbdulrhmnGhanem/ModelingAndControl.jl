% Chapter 5 Exercise 80
% Multisine response of a dynamic nonlinear system
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 6 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %% 
% Response nonlinear dynamic system to a periodic input %%
%                                                       %%
%   System considered                                   %%
%                                                       %%
%           y'' + 0.4*y' + y = g(u, u')                 %%
%   with                                                %%
%           g(u, u') = u*u'        (even nonlinearity)  %%
%           g(u, u') = u + u^2*u'  (odd nonlinearity)   %%
%                                                       %%
%   Inputs considered                                   %%
%       1. full multisine                               %%
%       2. odd multisine                                %%
%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition left hand side NL dynamic system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = 1;             % resonance angular frequency
zeta = 0.2;         % damping ratio


%%%%%%%%%%%%%%%%%%%%%%%%%
% multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%

% basic parameters multisine
f0 = 1;                                 % frequency multisine
fs = 200;                               % sampling frequency
N = fs/f0;                              % number of points in one period
t = [0:N-1].';                          % time in samples

% full multisine with 5 components: [1:5] * f0
Uf = zeros(N, 1) + eps;                 % add eps to make the non-excited lines visible in a dB-scale
kkf = [1:5].';
Uf(kkf+1) = exp(-sqrt(-1)*pi*kkf.*(kkf-1)/length(kkf));
uf = real(ifft(Uf));
uf = uf/std(uf, 1);                        % rms value equal to one

% odd multisine with 3 components: [1, 3, 5] * f0
Uo = zeros(N, 1) + eps;                 % add eps to make the non-excited lines visible in a dB-scale
kko = [1:2:5].';
Uo(kko+1) = exp(-sqrt(-1)*pi*kko.*(kko-1)/length(kko));
uo = real(ifft(Uo));
uo = uo/std(uo, 1);                            % rms value equal to one


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivative full multisine
Uf = fft(uf);
dUf = zeros(N, 1);
dUf(kkf+1) = Uf(kkf+1).*(2*pi*f0*sqrt(-1)*kkf);
duf = 2*real(ifft(dUf));

% derivative odd multisine
Uo = fft(uo);
dUo = zeros(N, 1);
dUo(kko+1) = Uo(kko+1).*(2*pi*f0*sqrt(-1)*kko);
duo = 2*real(ifft(dUo));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response to the nonlinear dynamic systems %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% odd input g(u, du/dt) = u + u^2 * du/dt second order differential equation
zOf = uf + uf.^2 .* duf;                % full multisine     
ZOf = fft(zOf);
zOo = uo + uo.^2 .* duo;                % odd multisine
ZOo = fft(zOo);

% even input  g(u, du/dt) = u * du/dt second order differential equation
zEf = uf .* duf;                        % full multisine     
ZEf = fft(zEf);
zEo = uo .* duo;                        % odd multisine
ZEo = fft(zEo);

% DFT lines at which the output is calculated (frequencies do not exceed 3*f0)
kk = [0:1:20].';

% output spectra
YOf = zeros(N,1);
YOo = zeros(N,1);
YEf = zeros(N,1);
YEo = zeros(N,1);

% transfer function
G = 1 ./ polyval([1, 2*zeta*w0, w0^2], 2*pi*f0*sqrt(-1)*kk);

% odd input nonlinearity, full multisine
YOf(kk+1) = G .* ZOf(kk+1);
yOf = 2*real(ifft(YOf));
% odd input nonlinearity, odd multisine
YOo(kk+1) = G .* ZOo(kk+1);
yOo = 2*real(ifft(YOo));
% even input nonlinearity, full multisine
YEf(kk+1) = G .* ZEf(kk+1);
yEf = 2*real(ifft(YEf));
% even input nonlinearity, odd multisine
YEo(kk+1) = G .* ZEo(kk+1);
yEo = 2*real(ifft(YEo));

% selection DFT lines shown in the figures
DFTlines = [0:20];
Uf = Uf(DFTlines+1)/N + eps;             % scaling by N to obtain the Fourier coefficients
Uo = Uo(DFTlines+1)/N + eps;             % add eps to show the non-excited lines in a dB scale
YOf = YOf(DFTlines+1)/N + eps;
YOo = YOo(DFTlines+1)/N + eps;
YEf = YEf(DFTlines+1)/N + eps;
YEo = YEo(DFTlines+1)/N + eps;


%%%%%%%%%%%%%%%%%%%%
% Plot the results %
%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(221)
plot(t, uf)
axis([-10, 210, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, yEf)
axis([-10, 210, -0.11, 0.11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uf), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(YEf), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Even Nonlinear Dynamic System - Full Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(2)
subplot(221)
plot(t, uo)
axis([-10, 210, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, yEo)
axis([-10, 210, -0.11, 0.11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uo), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(YEo), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Even Nonlinear Dynamic System - Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(3)
subplot(221)
plot(t, uf)
axis([-10, 210, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, yOf)
axis([-10, 210, -0.11, 0.11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uf), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(YOf), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Odd Nonlinear Dynamic System - Full Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

figure(4)
subplot(221)
plot(t, uo)
axis([-10, 210, -2.2, 2.2]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(t, yOo)
axis([-10, 210, -0.11, 0.11]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(DFTlines, db(Uo), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Input (dB)')

subplot(224)
plot(DFTlines, db(YOo), '+')
axis([-1, 21, -400, 200]);
xlabel('DFT line number')
ylabel('Output (dB)')

annotation('textbox',[0.2 0.8 0.05 0.2],'Color','k','String','Odd Nonlinear Dynamic System - Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

