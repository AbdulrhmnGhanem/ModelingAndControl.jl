% Chapter 1 Exercise 5a
% combining measurements with a varying SNR: Weighted Least squares
% estimation
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
NRepeat=100000              % Number of repeated experiments
NData=100                   % number of measurements per voltmeter in 1 experiment
iMax=0.01                   % maximum amplitude of the current
stdU1=1                     % std dev 1st voltmeter
stdU2=4                     % std dev 2nd voltmeter
R0=1000                     % exact value of the resistor

w=[ones(1,NData)*stdU1.^2 ones(1,NData)*stdU2.^2]';   % weight vector

for k=NRepeat:-1:1                      % loop over the repeated experiments
    i1=iMax*2*(rand(NData,1)-0.5);      % generate random current first experiment
    i2=iMax*2*(rand(NData,1)-0.5);      % generate random current 2nd experiment
    u1=R0*i1+stdU1*(randn(size(i1)));   % generate disturbed voltage
    u2=R0*i2+stdU2*(randn(size(i2)));   % generate disturbed voltage
    u=[u1' u2']';                       % combine data of both voltmeters
    i=[i1' i2']';
    R2(k)=(u'*(i./w))/(i'*(i./w));      % WLS
    R1(k)=(u'*i)/(i'*i);                %  LS
end


% calculate pdf
BinWidth=2;                             % set the parameters for the HIST call
X=[0:2:1200];

Rpdf(:,1)=hist(R1,X)/BinWidth/NRepeat;  % estimate pdf
Rpdf(:,2)=hist(R2,X)/BinWidth/NRepeat;  % estimate pdf

% plot the results
FigNum=6;

figure(FigNum),clf

plot(X,Rpdf(:,1),'k',X,Rpdf(:,2),'k')
axis([850 1150 0 0.03]);
DG_SetFontSize(7)
DG_SetTraceWidth(1.5,'*',FigNum)
DG_SetTraceTint(30,2,FigNum)
title('pdf least squares estimate')
xlabel('R')

% export the figure
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent4.pdf', gcf);

