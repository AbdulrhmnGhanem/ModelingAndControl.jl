% Chapter 7 Exercise 100
% Estimating a parametric model for Gbla using the robust method
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 8 December 2010

% in this routine, the a time domain object is prepared to be 
% loaded in the GUI of the FDIDENT-toolbox for further processing

clear all
N=1000;   % length of a period
M=2;      % Number of periods
R=7;      % Number of realizations of the input

fMaxGen=0.45;   % sets the bandwidth of the generator signal

NTrans=500;   % length transient part to eliminate transients in simulation


OrderG1=2 ;     % order of the plant G1
OrderG2=2;      % Order of the plant G2
NLpol=[1 1 1];  % coeff. of the static nonlinear polynomial function f
                % x, x^2, x^3
fMax1=0.25;     % to set the bandwidth of the first systems (fs=1)
fMax2=0.4;     % to set the bandwidth of the first systems (fs=1)
Amp=0.7;      % set the amplitude of the multisines

fPlot=[0:0.001:0.5];   % frequencies used to make the plots
Nnorm=floor(length(fPlot)*0.8); % used to normalize the gain of the estimated transfer function

% define the systems
[bGen,aGen]=butter(3,2*fMaxGen);bGen(end)=bGen(end);    % filter used to generate the excitation of the system
[b10,a10]=butter(OrderG1,2*fMax1); % plant to be estimated
[b20,a20]=cheby1(OrderG2,10,2*fMax1); % plant to be estimated

StdNoise=0.01;   % std. dev. of the disturbing noise


% Transfer functions for the plots
wPlot=fPlot*2*pi;   %freq grid
G10=freqz(b10,a10,wPlot);
G20=freqz(b20,a20,wPlot);
GBLA0=G10.*G20;    % GBLA
GGen=freqz(bGen,aGen,wPlot);

% create the data
for k=1:R    % loop over the realizations

% create the data
Lines=[1:1:floor(fMaxGen*N)];    % Lines to be excited
u0=RandMulti(Lines,N,N*M+NTrans)*Amp;   % generate random phase multisine as input

p0=filter(b10,a10,u0);                          % first linear system G1
q0=NLpol(1)*p0+NLpol(2)*p0.^2+NLpol(3)*p0.^3;   % static nonlinearity f
y0=filter(b20,a20,q0);                          % 2nd linear system G1

y=y0+StdNoise*randn(size(y0)); % add disturbing noise the output     

u0(1:NTrans)=[];y(1:NTrans)=[];   % eliminate transients
u0=u0(:);y=y(:);                % make collumn vector

uAll{1,k}=u0;yAll{1,k}=y;           % store individual realizations in a cell array


end

tdData=tiddata(yAll,uAll,1);
tdData.synchronization='samepower'
tdData.periodlength=N; 




