% Chapter 4 Exercise 58
% Identification in the time domain
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
na=3;                                   % order of the estimated model: denominator
nb=3;                                   %                               numerator
[b0,a0]=cheby1(3,5,2*0.1);              % system G0
[bGen,aGen]=cheby1(5,1,2*0.45);         % generation filter

% generate the signals
u0=filter(bGen,aGen,randn(N+NTrans,1));         % input signal
y=filter(b0,a0,u0)+StdNoise*randn(N+NTrans,1);  % output signal
u0(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate the simulation transients

% Part 1: set up the normal equations for ARX modelling, eliminating the initial state effects
n=max(na,nb);
yRHS=-y(n+1:end);                       % right handside of the ARX-normal equations
K=[toeplitz(y(n:end-1),y(n:-1:n-na+1))   -toeplitz(u0(n+1:end),u0(n+1:-1:n+1-nb))];  % K-matrix
Theta=K\yRHS;                           % least squares solution K Theta = yHRS
aEst=[1 Theta(1:na)'];                  % estimated model parameters: denominator
bEst=Theta(na+1:na+nb+1)';              %                             numerator

% Part 2: set up the normal equations for ARX modelling, without eliminating the initial state effects
u0Init=[zeros(n,1);u0];yInit=[zeros(n,1);y];    % extend with zeros
yRHS=-yInit(n+1:end);                           % right handside of the equations
KInit=[toeplitz(yInit(n:end-1),yInit(n:-1:n-na+1))  ...
          -toeplitz(u0Init(n+1:end),u0Init(n+1:-1:n+1-nb))];  % K-matrix
Theta=KInit\yRHS;                       % least squares solution K Theta = yHRS
aEstInit=[1 Theta(1:na)'];              % estimated model parameters: denominator
bEstInit=Theta(na+1:na+nb+1)';          %                             numerator

% plot the results
FigNum=1;
figure(FigNum)
clf

f=[0:500]/1000;                                 % frequencies for plot
[GEst]=freqz(bEst,aEst,2*pi*f);                 % estimated transfer function
[GEstInit]=freqz(bEstInit,aEstInit,2*pi*f);     % estimated transfer function
G0=freqz(b0,a0,2*pi*f);                         % exact transfer function

plot(f,db(GEst),'k',f,db(G0),'k',f,db(G0-GEst),'k',f,db(G0-GEstInit),'k')
DG_SetFontSize(7)
ylabel('Amplitude (dB)')
xlabel('frequency')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(1.5,1,FigNum,1)
DG_SetTraceWidth(1.5,3,FigNum,1)
DG_SetTraceWidth(0.5,2,FigNum,1)
DG_SetTraceTint(30,4,FigNum,1)
DG_SetTraceWidth(1.5,4,FigNum,1)
 
% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxTDnoNoiseNoModelErr.pdf', gcf);

