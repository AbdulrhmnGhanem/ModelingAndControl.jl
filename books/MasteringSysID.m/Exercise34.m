% Chapter 3 Exercise 34
% Impulse response function measurements
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

N=128               % Number of data points in the measurement
Ampl=1              % amplitude of the impulse
fs=256              % sample frequency

t=[0:N-1]/fs;       % time vector
f=[0:N-1]/N*fs;     % frequency vector
Lines=[1:N/2];      % select the lines to be plotted

[b,a]=cheby1(2,10,0.2);     % test system
G=freqz(b,a,f/fs*2*pi);     % calculate FRF

u=[Ampl zeros(1,N-1)]; % impulse to excite the system
y=filter(b,a,u);       % response of the system

U=fft(u)/sqrt(N);              % fft 
Y=fft(y)/sqrt(N);

% plot the results
FigNum=1
figure(FigNum)

DG_SetFontSize(12)

subplot(2,2,1)                  % impulse excitation
stem(t,u,'.k')
DG_SetFontSize(12)
xlabel('time (s)'),ylabel('u(t)')
title('Impulse excitation')
axis([-0.1 0.5 0 1.5])

subplot(2,2,2)                  % impulse response
stem(t,y,'.k')
DG_SetFontSize(12)
xlabel('time (s)'),ylabel('y(t)')
title('Impulse response')
axis([-0.1 0.5 -0.2 0.2])

subplot(2,2,3)                  % input spectrum
plot(f(Lines),abs(U(Lines)),'.k')
DG_SetFontSize(12)
title('Input spectrum'),xlabel('Frequency (Hz)'),ylabel('amplitude')
axis([0 0.5*fs 0 0.1])

subplot(2,2,4)                  % output spectrum
plot(f(Lines),abs(G(Lines)),'k',f(Lines),abs(Y(Lines)/U(1)),'.k')
DG_SetFontSize(12)
title('Output spectrum'),xlabel('Frequency (Hz)'),ylabel('amplitude')
axis([0 0.5*fs 0 1.5])
DG_SetTraceTint(50,1,FigNum,4)
DG_SetTraceWidth(2,1,FigNum,4)

% export the plots
%DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
%DG_MakePDF('FigExcitation301.pdf', gcf);    



