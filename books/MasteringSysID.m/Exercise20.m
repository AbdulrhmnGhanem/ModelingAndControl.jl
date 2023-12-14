% Chapter 2 Exercise 20
% Generation of a multisine with flat amplitude spectrum
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all
% define the input parameters
fsample=1000;Ts=1/fsample;      % sample frequency and sample period
Ndata1=1000;                    % length of the signal
t=[0:Ndata1-1]*Ts;              % time axis
Nsines=100;                     % number of sines
f=[0:Ndata1-1]*fsample/Ndata1;  % multisine frequencies
LinesPlot=[1:floor(Ndata1/2)];  % lines to be plotted

% multisine1: with a linear phase
U1=zeros(Ndata1,1);             % create spectrum 
tau=0.3;    
phiU1=-tau*2*pi*f;              % multisine with linear phase spectrum
U1(2:Nsines+1)=exp(j*phiU1(2:Nsines+1));     % choose all components in phase
u1=2*real(ifft(U1));u1=u1/std(u1);
U1m=fft(u1)/sqrt(Ndata1);       % spectrum of the actual generate multisine

% multisine 2: with a random phase
U2=zeros(Ndata1,1);             % choose random phases
U2(2:Nsines+1)=exp(j*2*pi*rand(Nsines,1));
u2=2*real(ifft(U2));u2=u2/std(u2);
U2m=fft(u2)/sqrt(Ndata1);       % spectrum of the actual generate multisine

% multisine 3: with a Schroeder phase
U3=zeros(Ndata1,1);
phse=-[1:Nsines].*[0:Nsines-1]*pi/Nsines;   % Schroeder phases
U3(2:Nsines+1)=exp(j*phse);
u3=2*real(ifft(U3));u3=u3/std(u3);
U3m=fft(u3)/sqrt(Ndata1);      % spectrum of the actual generate multisine

% plot the results
FigNum=5                        % time domain signals
figure(FigNum)
clf

subplot(1,2,1)
plot(t,u2,'k',t,u1,'k')
axis([0 1. -5 15])
DG_SetFontSize(9,FigNum)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,2,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0,FigNum)
xlabel('Time (s)')

subplot(1,2,2)
plot(t,u2,'k',t,u3,'k')
axis([0 1. -5 15])
DG_SetFontSize(9)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,2,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0,FigNum)
xlabel('Time (s)')

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigMultisine1a.pdf', gcf);

FigNum=FigNum+1;                % frequency domain plot
figure(FigNum)
plot(f(LinesPlot),db(U1m(LinesPlot)),'+k',...
     f(LinesPlot),db(U1m(LinesPlot)),'ok',...
     f(LinesPlot),db(U1m(LinesPlot)),'*k')
axis([0 500 -60 20])
DG_SetFontSize(12,FigNum)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,2,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0)
xlabel('Frequency (Hz)'),ylabel('Amplitude (dB)')

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigMultisine1b.pdf', gcf);

