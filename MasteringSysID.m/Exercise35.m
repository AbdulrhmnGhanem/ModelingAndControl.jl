% Chapter 3 Exercise 35
% Study of the sine response of a linear system: transients and steady
% state
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

N=128               % Number of data points in the measurement
fs=256              % sampling frequency
fExcit=16           % frequency of the sine wave
Nperiod=fs/fExcit;  % number of points/period 
Ampl=1              % amplitude of the sine


t=[0:N-1]/fs;       % time vector
f=[0:N-1]/N*fs;     % frequency vector
Lines=[1:N/2];      % select the lines to be plotted

[b,a]=cheby1(2,10,0.2);     % test system

u=cos(t*2*pi*fExcit);       % sine excitation of the system
y=filter(b,a,u);            % response of the system

U=fft(u)/N;                 % FFT analysis 
Y=fft(y)/N;                 % FFT analysis
Lines=[1:N/2];              % FFT lines to be plotted
f=[0:N-1]/N*fs;             % frequency axis

yEndPeriod=y(end-Nperiod+1:end);                    % take out the last period
M=N/Nperiod;                                        % number of periods
ySteadyState=kron(ones(1,M),yEndPeriod);            % repeat the last period fSine times
transient=y-ySteadyState;                           % calculate the transient


% Plot the results
FigNum=2;
figure(FigNum)
clf

subplot(1,2,1)                      % transient response
plot(t,y,'.k',t,transient,'.k')
xlabel('Time (s)'),ylabel('y(t) and transient')
title('Transient response')
axis([0 1/4 -1. 1.])
DG_SetTraceTint(50,1,FigNum,1)
DG_SetFontSize(12)

subplot(1,2,2)                      % FFT analysis of the output
plot(f(Lines),db(U(Lines)),'+k',f(Lines),db(Y(Lines)),'.k')
DG_SetTraceTint(50,2,FigNum,2)
axis([0 128 -60 0])
xlabel('f (Hz)'),ylabel('Amplitude (dB)')
title('FFT analysis of the output')
DG_SetFontSize(12)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation302.pdf', gcf);







