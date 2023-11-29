% Chapter 4 Exercise 65
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

% ARX estimation in the frequency domain; eliminate initial state effects
na=2;                                   % order of the estimated model with model errors
nb=2;
n=max(na,nb)+1;

U0=fft(u0);Y=fft(y);                    % fft

w=[0:N/2]/N*2*pi;w=exp(-sqrt(-1)*w(:)); % freq. vector on half of the unit circle
Lines=floor([0.04*N:1:0.08*N]);         % FFT lines to select the frequency band for fit
U0F=U0(Lines);YF=Y(Lines);wF=w(Lines);  % select the frequency band 

% set up the equations    frequency domain with initial conditions added
clear K
YRHS=[real(YF(:));imag(YF(:))];         % right handside of the equations
for k=nb:-1:0;                          % set up the normal equations
    K(:,na+k+1)=(wF.^k).*U0F;   
end
for k=na:-1:1
    K(:,k)=-(wF.^k).*YF;
end

for k=max(na,nb)-1:-1:0                 % add initial conditions
    K(:,na+nb+1+k+1)=(wF.^k);   
end

K=[real(K);imag(K)];
Theta=K\YRHS;                           % least squares solution K Theta = yHRS
aEstFD=[1 Theta(1:na)'];                % estimates: the denominator parameters
bEstFD=Theta(na+1:na+nb+1)';            %                nominator parameters

% plot the results
FigNum=10;
figure(FigNum)
clf

f=[0:N/2]/N;                             % frequency axis for plot
[GEstFD]=freqz(bEstFD,aEstFD,2*pi*f);    % FRF estimated transfer function
G0=freqz(b0,a0,2*pi*f);                  % FRF exact transfer function

plot(f,db(GEstFD),'k',f(Lines),db(GEstFD(Lines)),'k',f,db(G0),'k',...
             f,db(G0-GEstFD),'--k',f(Lines),db(G0(Lines)-GEstFD(Lines)),'--k')
axis([0 0.2 -40 5])
    
DG_SetFontSize(7)
ylabel('Amplitude (dB)')
xlabel('frequency')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceTint(30,5,FigNum,1)
DG_SetTraceWidth(1.5,'*',FigNum,1)
DG_SetTraceWidth(0.5,3,FigNum,1)

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxModErrBPFDG0Gest.pdf', gcf);  