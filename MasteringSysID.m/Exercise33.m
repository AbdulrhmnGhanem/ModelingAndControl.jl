% Chapter 2 Exercise 33
% Exploiting the periodic nature of signals: differentiation, integration,
% averaging and filtering
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters

N=1024                          % number of samples in a period
M=10                            % number of periods
SNR=0.001                       % scale the disturbing noise
fs=N                            % sampling period
[b,a]=butter(2,0.3);            % filter to reduce the bandwidth of the signal
t=[0:N-1]/fs;                   % time axis for plots

% generate M periods of a sawtooth signal
u0=[0:N/2-1 N/2:-1:1];u0=u0-mean(u0);u0=u0/std(u0); % 1 period normalized signal
u0=kron(ones(1,M+1),u0);                            % M+1 periodic repitition
u0=filter(b,a,u0);u0(1:N)=[];                       % Eliminate transient of the simulation
u=u0+SNR*randn(size(u0));                           % add white disturbing noise

% calculate spectrum and the derivative
n=length(u);
f=[0:n-1]/n*fs;                                     % frequency vector for plots
U=fft(u);U(n/2+1:end)=0;u1=2*real(ifft(U));         % original signal
dU=U.*f*sqrt(-1)*2*pi;  
du=2*real(ifft(dU));                                % derivative

U0=fft(u0);U0(n/2+1:end)=0;                         % fft noiseless signal 
dU0=U0.*f*sqrt(-1)*2*pi;
du0=2*real(ifft(dU0));du0=du0(1:N);                 % derivative of the noise les signal

dum=mean(reshape(du,N,M),2)';                       % average over the periods

dUm=fft(dum);dUm(1:2:end)=0;dumf=real(ifft(dUm));   % put equal frequency components to zero

FigNum=1
figure(FigNum)

subplot(2,2,1)
plot(t,du0,'k',t,u0(1:N),'k'),shg
axis([0 1 -20 20])
title('Noiseless (a)')
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,'*',FigNum,1)
DG_SetTraceWidth(1.5,1,FigNum,1)
DG_SetTraceTint(50,1,FigNum,1)

subplot(2,2,2)
plot(t,du(1:N),'k',t,du0,'k',t,u0(1:N),'k'),shg
axis([0 1 -20 20])
title('Noisy (b)')
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,'*',FigNum,2)
DG_SetTraceWidth(1.5,2,FigNum,2)
DG_SetTraceTint(50,2,FigNum,2)

subplot(2,2,3)
plot(t,dum,'k',t,du0,'k',t,u0(1:N),'k'),shg
axis([0 1 -20 20])
xlabel('Time (s)')
title('Averaged (c)')
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,'*',FigNum,3)
DG_SetTraceWidth(1.5,2,FigNum,3)
DG_SetTraceTint(50,2,FigNum,3)

subplot(2,2,4)
plot(t,dumf,'k',t,du0,'k',t,u0(1:N),'k'),shg
axis([0 1 -20 20])
xlabel('Time (s)')
title('Averaged and only odd freq. (d)')
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,'*',FigNum,4)
DG_SetTraceWidth(1.5,2,FigNum,4)
DG_SetTraceTint(50,2,FigNum,4)
DG_SetLineWidth(1.0)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Derivative1.pdf', gcf); 