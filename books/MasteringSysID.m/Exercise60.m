% Chapter 4 Exercise 60
% Numerical conditioning
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
N=5000                  % number of datapoints
NTrans=1000             % eliminate the transients of the simulation
StdNoise=0.             %  standard deviation of the disturbing noise put to zero

OrderAll=[1:15];                    % tested orders
BWAll=[1:2:68 70:1:100];            % scaling to reduce the bandwidth
BW0=[-0.2 0.2];                     % full bandwidth centered around 0.25  (fs=1)

r=randn(N+NTrans,1);                % reference signal unfiltered input

for k1=length(OrderAll):-1:1        % run over the orders
    for k2=length(BWAll):-1:1       % run over the bandwidth
        [k1 k2]                     % follow the for loops

    BW=0.25+BW0/BWAll(k2);          % select the bandwidth
    OrderG0=OrderAll(k1);           % select the order G0 system
    na=2*OrderG0;nb=2*OrderG0;      % order for the model, factor 2 for matlab code bandpass
        
    [b0,a0]=cheby1(OrderG0,20,2*BW);% system G0
    [bGen,aGen]=cheby1(5,5,2*BW);   % generation filter --> lower ripple, some power out of BW

% generate the signals
u0=filter(bGen,aGen,r);             % input signal
u0=u0/std(u0);                      % normalize
y=filter(b0,a0,u0)+StdNoise*randn(N+NTrans,1);  % output signal, no noise
u0(1:NTrans)=[];y(1:NTrans)=[];     % eliminate the transients of the simulation

% set up the normal equations for ARX modelling
n=max(na,nb)+1;

yRHS=-y(n+1:end);                   % right handside of the equations
K=[toeplitz(y(n:end-1),y(n:-1:n-na+1))   -toeplitz(u0(n+1:end),u0(n+1:-1:n+1-nb))];  % K-matrix
Theta=K\yRHS;                       % least squares solution K Theta = yHRS
aEst=[1 Theta(1:na)'];              % estimated model parameters
bEst=Theta(na+1:na+nb+1)';

CondAll(k1,k2)=ceil(log10(cond(K)));

% plot([db(freqz(b0,a0)) db(freqz(bEst,aEst))]),shg     % plot the intermediate results
  
    end   % loop over the bandwidth
end       % loop over the order


% plot the results
X=OrderAll(:);
Y=(BW0(2)-BW0(1))*2./BWAll(:);Y=log10(Y);

FigNum=1;
figure(FigNum),clf
[cs,h]=contour(Y,X,CondAll,[0:2:16],'k');clabel(cs,h,[1:1:16],'fontsize',6)
axis([-2 0 1 15]);shg
DG_SetFontSize(7,FigNum)
xlabel('log10(Relative Bandwidth)')
ylabel('Order')
DG_SetTraceWidth(1,'*',FigNum)  % Gm

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigARXCondNumber.pdf', gcf);  
