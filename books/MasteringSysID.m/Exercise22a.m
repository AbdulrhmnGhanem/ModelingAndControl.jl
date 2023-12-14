% Chapter 2 Exercise 22a
% Spectral analysis of a multisine, leakage present
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters
fSample=1000;                   % sample frequency
Ts=1/fSample;                   % sample period
Ndata1=128;                     % period length
M=4;                            % number of complete periods measured
Fraction=0.25;                  % fraction of the uncompletely measured period

% Generate the multisine
ExcitationLines=1+floor([0.15*Ndata1:0.35*Ndata1]);             % excited frequencies
U=zeros(Ndata1,1);                                              % initialize the spectrum
U(ExcitationLines)=exp(j*2*pi*rand(length(ExcitationLines),1)); % define the spectrum
u=2*real(ifft(U));                                              % generate 1 period
u=u/std(u);                                                     % normalize the signal to RMS=1
u=kron(ones(M,1),u);                                            % repeat M+1 times 
u(end-floor(Ndata1*(1-Fraction))+1:end)=[];                     % take out the last part to get the required length
Ndata=length(u);                                                % length of the signal
t=[0:Ndata-1]*Ts;                                               % time axis
f=[0:Ndata-1]*fSample/Ndata;                                    % frequency axis
LinesPlot=[1:Ndata/2];                                          % plot only half of the frequencies

% FFT analysis
URect=fft(u)/sqrt(Ndata);                                       % fft with rectangular window
HannWin=hanning(Ndata,'periodic');                              % fft with Hann window
UHann=fft(u.*HannWin(:))/sqrt(Ndata);

% plot the results
FigNum=7;
figure(FigNum)
clf

subplot(1,2,1)
plot(t,u,'k')
axis([0 0.6 -3 3])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0)
xlabel('Time (s)'),ylabel('Signal')

subplot(1,2,2)
plot(f(LinesPlot),db(URect(LinesPlot)),'.k',f(LinesPlot),db(UHann(LinesPlot)),'.k')
axis([0 500 -80 20])
xlabel('Frequency (Hz)')
ylabel('Amplitude Spectrum (dB)')
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,2,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0)

% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigMultisinus2.pdf', gcf);

