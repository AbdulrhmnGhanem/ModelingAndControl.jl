% Chapter 2 Exercise 25 and 26
% Generation of a Maximum length binary sequence
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameter
nRegister=5;                                        % register length
N=2^nRegister-1;                                    % period of the sequence
OverSample=16;                                      % oversampling to produce the plots
fSample=100;                                        % clock frequency

u=mlbs(nRegister);                                  % generate the maximum length sequence (FDiDent call)
uOverSample=kron(u(:),ones(OverSample,1));          % create oversampled signal

U=fft(u)/length(u);                                 % fft analysis
UOverSample=fft(uOverSample)/length(uOverSample);

f1=[0:N-1]'/N*fSample;                              % frequency vector for plots
f2=[0:N*OverSample-1]'/N*fSample;

% plot the results
FigNum=13
figure(FigNum)
clf

subplot(2,2,1)
stairs([0:N-1]'/fSample,u,'k')
axis([0 N/fSample -1.5 1.5])
hold on
plot([0:N-1]'/fSample,u,'.k')
hold off
xlabel('Time (s)')
title('ZOH-MLBS')
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,1,FigNum)
DG_SetLineWidth(1.0)

subplot(2,2,2)
stem(f1,abs(U(1:N)),'k')
xlabel('Frequency (Hz)')
ylabel('Amplitude (linear)')
title('FFT discrete time sequence')
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)

subplot(2,2,3)
stem(f1,abs(UOverSample(1:N)),'k');
xlabel('Frequency (Hz)')
ylabel('Amplitude (linear)')
title('FFT ZOH-MLBS')
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)

subplot(2,2,4)
plot(f2(1:N*OverSample/2),abs(UOverSample(1:N*OverSample/2)),'k')
xlabel('Frequency (Hz)')
ylabel('Amplitude (linear)')
title('FFT ZOH-MLBS')
DG_SetFontSize(12)%DG_SetLabelFormat('Y','%6.1f',FigNum,2)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('MLBS.pdf', gcf);
