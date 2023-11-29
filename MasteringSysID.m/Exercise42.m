% Chapter 3 Exercise 42
% Impulse response function measurements in the presence of output noise
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

N=128           % Number of data points in the measurement
NRepeat=100     % Number of repetitions to average
Ampl=1          % amplitude of the impulse
sny=0.02        % standard deviation noise in time domain
fs=256          % sample frequency
[b,a]=cheby1(2,10,0.2);  % test system


t=[0:N-1]/fs;   % time vector
f=[0:N-1]/N*fs; % frequency vector
Lines=[1:N/2];  % select the lines to be plotted

G=freqz(b,a,f/fs*2*pi);     % calculate FRF

u=[Ampl zeros(1,N-1)];      % impulse to excite the system
y0=filter(b,a,u(:));        % exact response of the system

for k=NRepeat:-1:1          % loop over the repeated experiments
    yNoise=sny*randn(N,1);
    y(:,k)=y0+yNoise;       % response of the system, disturbed with noise
end

U=fft(u)/N;                 % fft
Y0=fft(y0);
Y=fft(y);Ym=mean(Y,2);Ystd=std(Y,0,2)/sqrt(NRepeat);  % fft

% plot the results

FigNum=1
figure(FigNum)

DG_SetFontSize(14)
subplot(2,2,1)          % impulse excitations
stem(t,u,'.k')
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)
xlabel('Time (s)'),ylabel('u(t)')
title('Impulse excitation')
axis([-0.1 0.5 0 1.5])


subplot(2,2,2)          % impulse response
plot(t,yNoise,'k',t,y0(:,end),'k')
DG_SetFontSize(12)
DG_SetTraceTint(30,1,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
xlabel('Time (s)'),ylabel('y(t)')
title('Impulse response')
axis([-0.1 0.5 -0.2 0.2])

subplot(2,2,3)          % output spectrum single experiment 
plot(f(Lines),db(Y0(Lines,end)),'k',f(Lines),db(Y0(Lines)-Y(Lines,end)),'k+',...
                    f(Lines),db(Y(Lines)),'.k')   % plot last experiment
DG_SetFontSize(12)
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
title('Output spectrum single exp'),xlabel('Frequency (Hz)'),ylabel('Amplitude (dB)')
axis([0 0.5*fs -60 0.0])

subplot(2,2,4)          % output spectrum after averaging
plot(f(Lines),db(Y0(Lines,end)),'k',f(Lines),db(Ym(Lines)-Y0(Lines,end)),'k+',...
                    f(Lines),db(Ym(Lines)),'.k',f(Lines),db(Ystd(Lines)),'k')   % plot last experiment
DG_SetFontSize(12)
DG_SetTraceTint(30,2,FigNum,4)
DG_SetTraceTint(60,4,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
title('Averaged output spectrum'),xlabel('Frequency (Hz)'),ylabel('Amplitude (dB)')
axis([0 0.5*fs -60 0])

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation311.pdf', gcf);  



