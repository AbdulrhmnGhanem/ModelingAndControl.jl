% Chapter 2 Exercise 16
% Generate a sine wave, no integer number of periods measured
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input variables
Ndata1=16                       % Number of data points in the first sine wave
Ndata2=2*Ndata1                 % Number of data points in the second sine wave
Ampl = 1                        % amplitude of the sine
phse=pi/2                       % phase of the sine
fSample=1000                    % sample frequency
fSine=80                        % frequency of the sinewave
Ts=1/fSample;                   % sample period
t=[0:Ndata1-1]*Ts;              % time vector


u=Ampl*sin(2*pi*fSine*t+phse);  % time signal
U=fft(u)/length(u);             % DFT

Lines1=[1:Ndata1]-1;            % DFT line numbers
f1=(Lines1)*fSample/Ndata1;     % frequency axis

% plot the results
FigNum=1
figure(FigNum)

subplot(2,2,1)
plot(t,u,'.k')
                            
subplot(2,2,2)
plot(f1,db(U),'+k')
axis([0 1000 -40 0])

subplot(2,2,3)
stem(Lines1,abs(U),'+k')
axis([0 16 0 1])
                             
subplot(2,2,4)
stem(f1,abs(U),'+k')
axis([0 1000 0 1])
                             
subplot(2,2,1)
DG_SetFontSize(12,FigNum,1)
DG_SetTraceWidth(0.5,1,FigNum,1)
xlabel('Time (s)'),ylabel('y(t)')

subplot(2,2,2)
DG_SetFontSize(12,FigNum,2)
DG_SetTraceWidth(0.5,1,FigNum,2)
title('DFT'),xlabel('Frequency (Hz)'),ylabel('Amplitude (dB)')

subplot(2,2,3)
DG_SetFontSize(12,FigNum,3)
DG_SetTraceWidth(0.5,1,FigNum,3)
title('DFT'),xlabel('DFT Line Number'),ylabel('Amplitude (Linear)')

subplot(2,2,4)
DG_SetFontSize(12,FigNum,4)
DG_SetTraceWidth(0.5,1,FigNum,4)
title('DFT'),xlabel('Frequency (Hz)'),ylabel('Amplitude (Linear)')

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigSinus1.pdf', gcf);

 