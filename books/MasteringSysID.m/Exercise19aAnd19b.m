% Chapter 2 Exercise 19.a and 19.b
% Generate a sine wave, using the Matlab IFFT instruction
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all
% define the input variables
Ndata1=16                       % Number of data points in the first sine wave
Ampl = 1                        % amplitude of the sine
phse=pi/2                       % phase of the sine
fSample=1000                    % sample frequency
fSine=fSample/Ndata1            % frequency of the sinewave, chosen to avoid leakage
Ts=1/fSample;                   % sample period
t=[0:Ndata1-1]*Ts;              % time vector


% Part 1: using the full spectrum
U1=zeros(Ndata1,1);             % create spectrum
U1(2)=exp(-j*phse);U1(Ndata1)=exp(+j*phse);  

u1=ifft(U1)*Ndata1/2;           % IFFT inverse Fourier transform

% Part 2: using only the positive frequencies
U2=zeros(Ndata1,1);    % create spectrum
U2(2)=exp(-j*pi/2); 

u2=Ampl*2*real(ifft(U2))*Ndata1/2;            % inverse Fourier transform

% plot the results
FigNum=4
figure(FigNum)
clf
plot(t,u1,'ok',t,u2,'+k')
DG_SetFontSize(9)
DG_SetTraceWidth(0.5,2,FigNum)
DG_SetTraceWidth(1,1,FigNum)
xlabel('Time (s)')
 
% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigSinus4.pdf', gcf);

