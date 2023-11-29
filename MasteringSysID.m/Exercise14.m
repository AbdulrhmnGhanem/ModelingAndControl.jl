% Chapter 2 Exercise 14
% Discretisation in time: choice of the sampling frequency: ALIAS
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 30 November 2010

clear all

% signal parameters: define a sine wave signals
A=1;                        % amplitude of the sine waves
f1=1;                       % frequency of the first sine in Hz
f2=9;                       % frequency of the 2nd sine in Hz

% measurement parameters
T=1;                        % measurement time in s
fs=10;                      % sample frequency in Hz

tc=[0:T/1000:T];            % time vector for continuous time plots
td=[0:1/fs:T-1/fs];         % time vector for discrete time plots
N=length(td);               % number of samples
fd=[0:N-1]*fs/N;            % frequencies of discrete Fourier transform

% create the signals
uc1=A*sin(2*pi*f1*tc);      % first 'continuous time' signal
uc2=-A*sin(2*pi*f2*tc);     % 2nd 

ud1=A*sin(2*pi*f1*td);      % first 'discrete time' signal
ud2=-A*sin(2*pi*f2*td);     % 2nd 

U1=fft(ud1)/N;              % DFT result uc1
U2=fft(ud2)/N;              % DFT result uc2

% plot the signals
FigNum=1;
figure(FigNum),clf

subplot(3,2,1)   % plot slow sine + samples
plot(tc,uc1,'-k',td,ud1,'ok')
ylabel('Slow sine')                              
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,1)   
DG_SetTraceTint(50,1,FigNum,1)

subplot(3,2,2)     % plot spectrum slow sine
stem([-f1 f1],[A/2 A/2],'filled','k')
hold on,stem([-fs fs],[0.93 0.93],'^:k'),hold off
axis([-15 15 0 1])
ylabel('Amplitude')
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,2)   
DG_SetTraceTint(50,2,FigNum,2)

subplot(3,2,3)     % plot sampled slow and fast sine
plot(td,ud1,'ok',td,ud2,'*k')
ylabel('Sampled sines')
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,3)   

subplot(3,2,4)   % plot dft spectrum
stem(fd,abs(U1),'ok'),hold on,stem(fd,abs(U2),'*k'),hold off
axis([-15 15 0 1])
ylabel('Amplitude')
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,4)   


subplot(3,2,5)   % plot continuous time sine + samples
plot(tc,uc2,'-k',td,ud2,'*k')
xlabel('Time (s)')
ylabel('Fast sine')
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,5)   
DG_SetTraceTint(50,1,FigNum,5)



subplot(3,2,6)     % plot spectrum continuous time fast sine
stem([-f2 f2],[A/2 A/2],'filled','k')
hold on,stem([-fs fs],[0.93 0.93],'^:k'),hold off
axis([-15 15 0 1])
ylabel('Amplitude')
xlabel('Frequency (Hz)')
DG_SetFontSize(12)
DG_SetLabelFormat('Y','%6.1f')
DG_SetTraceWidth(0.5,1,FigNum,6)   


% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigDiscretAlias.pdf', gcf);
%print -dpdf FigDiscretAlias