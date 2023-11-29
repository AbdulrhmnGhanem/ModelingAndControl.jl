% Chapter 2 Exercise 29
% Repeated realizations of a white random noise excitation with increasing
% length
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

% define the input parameters
Ndata=128                           % Number of data points in the random excitation
Ampl=1                              % RMS value of the signal
fSample=1000                        % sample frequency

% Generate 4 realizations with the same length
FigNum=2
figure(FigNum)

for k=1:4
    N=Ndata*2^(k-1);                % change the length of the record
    Lines=1:N/2 ;                   % plotted FFT lines
    f=(Lines-1)/N*fSample;
    u=randn(N,1);                   % generate the random signal
    u=Ampl*u/std(u);                % set the rms value
    U=fft(u)/sqrt(N);               % fft analysis
   
    subplot(2,2,k),plot(f,db(U(Lines)),'.k')
    axis([0 500 -30 10])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    title(['Record length ',int2str(N)])
    DG_SetFontSize(10)
    DG_SetTraceWidth(1,1,FigNum)
    DG_SetLineWidth(1.0)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Noise2.pdf', gcf); 