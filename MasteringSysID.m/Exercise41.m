% Chapter 3 Exercise 41
% FRF-measurements using a burst excitation
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables
N=256                           % Number of data points in the measurement
NBurstAll=256-[64*3 64*2 64];   % different burst lengths
fs=N                            % sampling frequency s.t. 1 period = 1 s
Ampl=1                          % RMS-value of the burst

[b,a]=cheby1(2,20,0.2);         % test system with high resonance
b=b*10;                         %                increased gain

t=[0:N-1]/fs;                   % time vector 
Lines=[1:N/2];                  % select the lines to be plotted
f=[Lines-1]'/N*fs;              % frequency vector
G0=freqz(b,a,f*2*pi/fs);        % exact FRF




for k=length(NBurstAll):-1:1       % loop over the burst lengths
       
% generate excitation signal
    u(:,k)=zeros(N,1);
    u(1:NBurstAll(k),k)=randn(NBurstAll(k),1)*Ampl;   % white noise excitation
    
% system response    
    y(:,k)=filter(b,a,u(:,k));                 % response of the system

Y=fft(y(:,k))/sqrt(N);           % spectral analysis
U=fft(u(:,k))/sqrt(N);

G(:,k)=Y(Lines)./U(Lines);                     % FRF

end


% plot the reuslts 
%        upper row: the time signals
%        lower row: the FRF and the errors

FigNum=9;
figure(FigNum)
clf

subplot(2,3,1)
plot(t,u(:,1),'k',t,y(:,1),'k')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,3,4)
plot(f,db(G(:,1)),'k',f,db(G0),'k',f,db(G(:,1)-G0),'.k')
ylabel('Amplitude (dB)')
axis([0 60 -40 20])
DG_SetTraceTint(30,1,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,3,2)
plot(t,u(:,2),'k',t,y(:,2),'k')
xlabel('Time (s)')
title('Input/Output')
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,3,5)
plot(f,db(G(:,2)),'k',f,db(G0),'k',f,db(G(:,2)-G0),'.k')
axis([0 60 -40 20])
title('FRF')
DG_SetTraceTint(30,1,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
xlabel('Frequency (Hz)')
DG_SetFontSize(12)

subplot(2,3,3)
plot(t,u(:,3),'k',t,y(:,3),'k')
DG_SetTraceTint(30,2,FigNum,5)
DG_SetTraceWidth(2,'*',FigNum,5)
DG_SetFontSize(12)

subplot(2,3,6)
plot(f,db(G(:,3)),'k',f,db(G0),'k',f,db(G(:,3)-G0),'.k')
axis([0 60 -40 20])

DG_SetTraceTint(30,1,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)

DG_SetFontSize(12)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation309.pdf', gcf);  
  
  
  