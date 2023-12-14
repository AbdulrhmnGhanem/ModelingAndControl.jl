% Chapter 4 Exercise 64
% Shaping the model errors in the time domain: pre-filtering
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
Lines=floor([0.04*N:1:0.08*N]);         % DFT lines used for plotting

[b0,a0]=cheby1(5,5,2*0.1);              % system G0
[bGen,aGen]=cheby1(5,1,2*0.15*3);       % generation filter

% generate the signals
u0=filter(bGen,aGen,randn(N+NTrans,1));         % input signal
y=filter(b0,a0,u0)+StdNoise*randn(N+NTrans,1);  % output signal
u0(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate the simulation transients

[bf,af]=butter(5,[0.04,0.08]*2);        % create prefilter
u0f=filter(bf,af,u0);                   % prefilter the input
yf=filter(bf,af,y);                     %           the output

% ARX estimation in the time domain; eliminate initial state effects
na=2;                                   % order of the estimated model with model errors
nb=2;
n=max(na,nb)+1;

yRHS=-yf(n+1:end);                      % right handside of the equations
K=[toeplitz(yf(n:end-1),yf(n:-1:n-na+1))   -toeplitz(u0f(n+1:end),u0f(n+1:-1:n+1-nb))];  % K-matrix
Theta=K\yRHS;                           % least squares solution K Theta = yHRS
aEst=[1 Theta(1:na)'];                  % estimated model parameters: denominator
bEst=Theta(na+1:na+nb+1)';              %                             numerator


% plot the results
FigNum=9;
figure(FigNum)
clf

f=[0:N/2]/N;                            % frequencies for plot
[GEst]=freqz(bEst,aEst,2*pi*f);         % estimated transfer function
G0=freqz(b0,a0,2*pi*f);                 % exact transfer function
fpref=[0:0.001:0.2];Gpref=freqz(bf,af,2*pi*fpref);  % G prefilter

plot(f,db(GEst),'k',f(Lines),db(GEst(Lines)),'k',f,db(G0),'k',...
             f,db(G0-GEst),'--k',f(Lines),db(G0(Lines)-GEst(Lines)),'--k',...
             fpref,db(Gpref),'k')
axis([0 0.2 -40 5])
    
DG_SetFontSize(7)
ylabel('Amplitude (dB)')
xlabel('frequency')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceTint(30,5,FigNum,1)
DG_SetTraceTint(20,6,FigNum,1)
DG_SetTraceWidth(1.5,'*',FigNum,1)
DG_SetTraceWidth(0.5,3,FigNum,1)
DG_SetTraceWidth(0.75,6,FigNum,1)

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxModErrBPTDG0Gest.pdf', gcf);  
  