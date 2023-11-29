% Chapter 2 Exercise 21
% The swept sine signal
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all
% define the input parameters
fsample=1000;                       % sample frequency
Ts=1/fsample;                       % sample period
NPeriod=1000;                       % length of a period
t=[0:NPeriod-1]*Ts;                 % time axis
k1=50;k2=200;                       % define the frequency band to be excited
f=[0:NPeriod/2-1]*fsample/NPeriod;  % swept sine frequencies for plot

% swept sine parameters
f0=1/(NPeriod*Ts);
a=pi*(k2-k1)*f0^2;b=2*pi*k1*f0;

u=sin((a*t+b).*t);                   % generate the swept sine

U=fft(u(1:NPeriod))/sqrt(NPeriod);   % calculate the FFT

% plot the results
FigNum=12
figure(FigNum)
clf

subplot(1,2,1)
plot(t,u,'k')
axis([0 1 -1.1 1.1])
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0,FigNum)
xlabel('Time (s)')

subplot(1,2,2)
plot(f,db(U(1:NPeriod/2)),'k')
axis([0 500 -50 10])
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0,FigNum)
xlabel('Frequency (Hz)')

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('SweptSine.pdf', gcf);
