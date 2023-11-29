% Chapter 2 Exercise 15
% Windowing: study of the leakage effect and the frequency resolution
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 30 November 2010

% signal parameters: define a sine wave signals
A=1;                        % amplitude of the sine waves
f1=2;                       % frequency of the first sine in Hz

% measurement parameters
T=3;                        % visualisation time in s of continous time signal
fs=10;                      % sample frequency in Hz
NAll=[10 11 12 13 14 15];   % number of samples in discrete time signal

tc=[0:T/1000:T];            % time vector for continuous time plots
uc=A*sin(2*pi*f1*tc);       % 'continuous time' signal

FigNum=1;
figure(FigNum),clf
n=length(NAll);

for k=1:n                   % make plots for different data lengths
    N=NAll(k);              % actual number of points
    td=0.5+[0:N-1]/fs;      % time vector for discrete time plots
    tw=[td(1):1/(fs*100):td(end)+1/fs];  % actual time window
    uw=A*sin(2*pi*f1*tw);   % continuous time signal in the window
    
    fd=[0:N-1]*fs/N;        % frequencies of discrete Fourier transform
    ud=A*sin(2*pi*f1*td);   % 'discrete time' signal

    U=fft(ud)/N;            % DFT result ud

    subplot(n,2,(k-1)*2+1)  % plot continuous sine + samples
    plot(tw,uw,'-k',tc,uc,'-k',td,ud,'ok')
    axis([0 T -1 1])
    
    DG_SetFontSize(12)
	DG_SetLabelFormat('Y','%6.1f')
	DG_SetTraceWidth(0.5,2,FigNum,(k-1)*2+1)
    DG_SetTraceWidth(2,1,FigNum,(k-1)*2+1)
    DG_SetTraceTint(50,2,FigNum,(k-1)*2+1)
    DG_SetLabelFormat('Y','%6.1f')

	subplot(n,2,(k-1)*2+2)     % plot spectrum slow sine
	stem([f1],[ A/2],'filled','k')
	hold on,stem(fd,abs(U),'+:k'),hold off
	axis([0 10 0 1])

	DG_SetFontSize(12)
	DG_SetLabelFormat('Y','%6.1f')
	DG_SetTraceWidth(0.5,1,FigNum,(k-1)*2+2)
    DG_SetTraceWidth(2,3,FigNum,(k-1)*2+2)  
end

subplot(n,2,(k-1)*2+1)
xlabel('Time (s)')

subplot(n,2,(k-1)*2+2)   
xlabel('Frequency (Hz)')

% DG_Init4PDF(gcf, 4.5, 4.5/1.618/2);        % fixing the size, half standard height
% DG_MakePDF('FigWindowLeakage.pdf', gcf);

