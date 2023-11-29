% Chapter 4 Exercise 76
% Using the frequency domain identification toolbox FDIDENT
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
NPer=1024;                              %  period length
M=7;                                    % number of processed periods
MTrans=1;MAll=M+MTrans;                 % number of periods for transient elimination of the simulation
NTrans=NPer*MTrans;                     % transient points
N=M*NPer;                               % total number of processed data
NRep=10;                                % number of repeated simulations
stdNoise=0.1;                           % standard deviation of disturbing noise
OrderG=2;                               % order of the plant
OrderNoise=1;                           % order of the noise filter

% plant model
[b0,a0]=cheby1(OrderG,10,2*0.25);b0(2)=b0(2)*1.3; % plant to be estimated

% noise model
[b,a]=butter(OrderNoise,2*0.2); b=b+0.001*a;% noise model
bNoise=b;   
aNoise=a; 

% generate the signals
Lines=[1:NPer/3];
u0=RandMulti(Lines,NPer,NPer*MAll);     % flat random phase multisine excitation
y0=filter(b0,a0,u0);                    % exact output
yNoise=filter(bNoise,aNoise,randn(1,N+NTrans));
yNoise=stdNoise*yNoise;                 % disturbing noise
y=y0+yNoise;                            % disturbed output

u0(1:NPer)=[];y(1:NPer)=[];  % eliminate transients
% create an object
ExpData=tiddata(y(:),u0(:),1);
   
   