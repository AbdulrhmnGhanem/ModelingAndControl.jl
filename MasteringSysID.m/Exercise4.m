% Chapter 1 Exercise 4
% Importance of the choice of the independent variable or input
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010
%

clear all
NRepeat=100000              % Number of repeated experiments
NData=100                   % number of measurements in 1 experiment
iMax=0.01                   % maximum amplitude of the current
stdU=1                      % standard deviation noise on the voltage
stdIAll=0                   % no noise on the current
R0=1000                     % exact value of the resistor

for r=length(stdIAll):-1:1  % loop over the std. dev. of the current
    stdI=stdIAll(r);
	for k=NRepeat:-1:1      % loop over the repeated experiments
        i=iMax*2*(rand(NData,1)-0.5);  % generate exact random current excitation 
        uN=R0*i+stdU*randn(size(i));   % generate normaly disturbed output
        iN=i+stdI*randn(size(i));      % add noise to the current
        R1(k,r)=uN'*i/(i'*i);          % estimate R with i as regressor
        R2(k,r)=1./(uN'*iN/(uN'*uN));  % estimate R with u as regressor
    end
end

% calculate pdf
BinWidth=2;     % set the parameters for HIST call
X=[0:2:1100];

Rpdf(:,1)=hist(R1(:,1),X)/BinWidth/NRepeat;  % estimate pdf
Rpdf(:,2)=hist(R2(:,1),X)/BinWidth/NRepeat;  % estimate pdf


% plot the results
FigNum=5;

figure(FigNum)
clf

plot(X,Rpdf(:,1),'k',X,Rpdf(:,2),'k')
axis([900 1100 0 0.025])
DG_SetFontSize(7)
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(1.5,'*',FigNum,2*(k-1)+1)
    
title('pdf of the estimates')
xlabel('R')


% Export the results
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent3Tris.pdf', gcf);

