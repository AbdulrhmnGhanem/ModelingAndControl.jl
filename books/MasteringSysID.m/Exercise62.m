%  Identification of linear dynamic systems 
%          using an ARX model and no disturbing noise
%
% Chapter 4 Exercise 62
% Identify a too simple model
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010
clear all

% define the input variables
N=5000                                  % number of datapoints
NTrans=1000                             % eliminate the transients of the simulation
StdNoise=0.                             %  standard deviation of the disturbing noise put to zero

[b0,a0]=cheby1(5,5,2*0.1);              % system G0
[bGen,aGen]=cheby1(5,1,2*0.15*3);       % generation filter

% generate the signals
u0=filter(bGen,aGen,randn(N+NTrans,1));         % input signal
y=filter(b0,a0,u0)+StdNoise*randn(N+NTrans,1);  % output signal
u0(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate the simulation transients

% estimate the model
na=4;                                   % order of the estimated model with model errors
nb=4;
n=max(na,nb)+1;

yRHS=-y(n+1:end);                       % right handside of the normal equations
K=[toeplitz(y(n:end-1),y(n:-1:n-na+1))   -toeplitz(u0(n+1:end),u0(n+1:-1:n+1-nb))];  % K-matrix
Theta=K\yRHS;                           % least squares solution K Theta = yHRS
aEst=[1 Theta(1:na)'];                  % estimated model parameters: the denominator
bEst=Theta(na+1:na+nb+1)';              %                             the numerator

% plot the results
FigNum=7;
figure(FigNum)
clf

f=[0:500]/1000;                         % frequencies for the plot
[GEst]=freqz(bEst,aEst,2*pi*f);         % estimated transfer function
G0=freqz(b0,a0,2*pi*f);                 % exact transfer function

plot(f,db(GEst),'k',f,db(G0),'k')
axis([0 0.5 -80 20])
DG_SetFontSize(7)
ylabel('Amplitude (dB)')
xlabel('frequency')
DG_SetTraceTint(100,1,FigNum,1)
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(1.5,'*',FigNum,1)
  
% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxModErrTDG0Gest.pdf', gcf);  
