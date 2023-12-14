% Chapter 2 Exercise 28
% Repeated realizations of a white random noise excitation with fixed length
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters
Ndata=128                       % Number of data points in the random excitation
Ampl=1                          % RMS value of the signal
fSample=1000                    % sample frequency
Lines=1:Ndata/2 ;               % plotted FFT lines
f=(Lines-1)/Ndata*fSample;      % frequency vector

% Generate 4 realizations with the same length and plot

FigNum=1;
figure(FigNum)
clf

for k=1:4
    u=randn(Ndata,1);           % generate the random signal
    u=Ampl*u/std(u);            % set the rms value
    U=fft(u)/sqrt(Ndata);       % fft analysis
    
    subplot(2,2,k),plot(f,db(U(Lines)),'.k')
    axis([0 500 -30 10])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')   
    DG_SetFontSize(12)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetLineWidth(1.0)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Noise1.pdf', gcf);    
