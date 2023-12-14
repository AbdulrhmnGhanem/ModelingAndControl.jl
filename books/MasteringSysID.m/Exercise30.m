% Chapter 2 Exercise 30
% Smoothing the amplitude spectrum of a random excitation
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all
% define the input parameters

Ndata=2048                      % Number of data points in the random excitation
AmpBin=[0:sqrt(Ndata)]/sqrt(Ndata)*4; 
                                % bin centers to be used in histogram plot
Naverage=[1 4 16];              % Number of averages used
Ampl=1                          % RMS value of the signal
fSample=1000                    % sample frequency

Lines=1:Ndata/2 ;               % plotted FFT lines
f=(Lines-1)/Ndata*fSample;      % plotted frequencies
n=length(Naverage);             % number of averaged realizations

% Generate the signals and plot the results
FigNum=3;
figure(FigNum)
clf

for r=1:n
    Naver=Naverage(r);          % number of averages
    u=randn(Ndata,Naver);       % generate the random signals (1 collumn is 1 realization)
    U=fft(u)/sqrt(Ndata);       % fft analysis
    
    U=sqrt(mean(abs(U(Lines,:).^2),2));  % average over the different realizations
    
    subplot(2,n,r),plot(f,db(U),'.k')
    axis([0 500 -30 10])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    title(['Number of averages = ',int2str(Naver)])
    DG_SetFontSize(10)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetLineWidth(1.0)
    
    subplot(2,n,n+r)
    H=hist(U,AmpBin);
    bar(AmpBin,H,'k')
    xlabel('Amplitude (linear)')
    title('Histogram')
    axis([0 4 0 300])
    DG_SetFontSize(10)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetLineWidth(1.0)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Noise3.pdf', gcf); 


