% Chapter 4 Exercise 59
% Identification in the frequency domain
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

% Part 1: frequency domain implementation,  no elimination of initial state effects

U0=fft(u0);Y=fft(y);                    % fft
U0(N/2+2:end)=[];Y(N/2+2:end)=[];       % select positive frequency axis
w=[0:N/2]/N*2*pi;w=exp(-sqrt(-1)*w(:)); % freq. vector on half of the unit circle to set up normal equations

YRHS=[real(Y(:));imag(Y(:))];           % set up of the normal equations
for k=nb:-1:0;    
    K(:,na+k+1)=(w.^k).*U0;   
end
for k=na:-1:1
    K(:,k)=-(w.^k).*Y;
end
KReal=[real(K);imag(K)];                % split in real and imaginary part --> real least squares solution
Theta=KReal\YRHS;                       % least squares solution K Theta = yHRS
aEstFDnoInit=[1 Theta(1:na)'];          % estimated model parameters: denominator 
bEstFDnoInit=Theta(na+1:na+nb+1)';      %                             numerator

% Part 2: frequency domain implementation, elimination of the initial state
% effects

for k=max(na,nb)-1:-1:0                 % add initial state effects to the normal equations
    K(:,na+nb+1+k+1)=(w.^k);   
end

K=[real(K);imag(K)];                    % split in real and imaginary part --> real least squares solution
Theta=K\YRHS;                           % least squares solution K Theta = yHRS
aEstFDWithInit=[1 Theta(1:na)'];        % estimated model parameters: denominator 
bEstFDWithInit=Theta(na+1:na+nb+1)';    %                             numerator


% plot the results
FigNum=2;
figure(FigNum)
clf

f=[0:500]/1000;                         % frequencies for the plot
[GEstFDnoInit]=freqz(bEstFDnoInit,aEstFDnoInit,2*pi*f);             % estimated transfer function
[GEstFDWithInit]=freqz(bEstFDWithInit,aEstFDWithInit,2*pi*f);       % estimated transfer function
G0=freqz(b0,a0,2*pi*f);                                             % exact transfer function

plot(f,db(GEstFDWithInit),'k',f,db(G0),'k',f,db(G0-GEstFDWithInit),'k',f,db(G0-GEstFDnoInit),'k')

DG_SetFontSize(7)
ylabel('Amplitude (dB)')
xlabel('frequency')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceTint(30,4,FigNum,1)
DG_SetTraceWidth(1.5,1,FigNum,1)
DG_SetTraceWidth(1.5,3,FigNum,1)
DG_SetTraceWidth(1.5,4,FigNum,1)
DG_SetTraceWidth(0.5,2,FigNum,1)
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxFDnoNoiseNoModelErr.pdf', gcf);  

