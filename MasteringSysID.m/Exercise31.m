% Chapter 2 Exercise 31
% Generation of random noise excitations with a user imposed power spectrum
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters
Ndata=2048                          % Number of data points in the random excitation
Naverage=[1 16 64];                 % Number of averaged realizations
Ampl=1                              % RMS value of the signal
fSample=1000                        % sample frequency
CutOff=0.1                          % Cut off frequency of the filter, relatively to fSample
[b,a]=butter(6,CutOff*2);           % noise shaping filter

Lines=1:Ndata/2 ;                   % plotted FFT lines
f=(Lines-1)/Ndata*fSample;          % plotted frequencies
n=length(Naverage);                 % number of experiments

G0=db(freqz(b,a,2*pi*f/fSample));   % FRF of the desired filter

% Generate and analyse the signals
FigNum=4;
figure(FigNum)
clf

for r=1:n 
    Naver=Naverage(r);              % number of averages
    clear URect UHann
    for k=Naver:-1:1                % run over the averaging loop
        x=randn(Ndata,1);           % generate the filtered random signals 
        x=filter(b,a,x);            

        URect(:,k)=fft(x)/sqrt(Ndata);   % fft analysis with rectangular window
                                    % (1 collumn is 1 realization)
        UHann(:,k)=fft(hanning(Ndata,'periodic').*x)/sqrt(Ndata);
                                    % fft analysis with hanning window
    end
    
    URectm=sqrt(mean(abs(URect(Lines,:).^2),2));        % average over the different realizations
    Scale=sqrt(mean(hanning(Ndata,'periodic').^2));     % scale factor for the Hanning window
    UHannm=sqrt(mean(abs(UHann(Lines,:).^2),2))/Scale;
    
    % Plot the results
    subplot(2,n,r),plot(f,db(URectm),'.k',f,G0,'k')     % plot rectangular analysis
    axis([0 500 -90 10])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    title(['Number of averages = ',int2str(Naver)])
    DG_SetFontSize(12)
    DG_SetTraceWidth(0.5,2,FigNum)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetTraceGray(0.5,1,FigNum)
    DG_SetLineWidth(1.0)
    
    subplot(2,n,r+n),plot(f,db(UHannm),'.k',f,G0,'k')     % plot hanning analysis
    axis([0 500 -90 10])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    title(['Number of averages = ',int2str(Naver)])
    DG_SetFontSize(12)
    DG_SetTraceWidth(0.5,2,FigNum)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetTraceGray(0.5,1,FigNum)
    DG_SetLineWidth(1.0)

end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Noise4.pdf', gcf);     



