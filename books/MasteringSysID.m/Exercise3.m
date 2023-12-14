% Chapter 1 Exercise 3
% Impact of disturbing noise on the regressor or input measurements
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
NRepeat=100000              % Number of repeated experiments
NData=100                   % number of measurements in 1 experiment
iMax=0.01                   % maximum amplitude of the current
stdU=1                      % standard deviation noise on the voltage
stdIAll=iMax*[0:0.5:1]/10   % standard deviation noise on the current
R0=1000                     % exact value of the resistor

for r=length(stdIAll):-1:1  % loop over the std. dev. of the current
    stdI=stdIAll(r);
	for k=NRepeat:-1:1      % loop over the repeated experiments
        i=iMax*2*(rand(NData,1)-0.5);  % generate exact random current excitation
        uN=R0*i+stdU*randn(size(i));   % generate normaly disturbed output
        iN=i+stdI*randn(size(i));      % add noise to the current
        R1(k,r)=uN'*i/(i'*i);          % estimate R, no noise on current
        R2(k,r)=uN'*iN/(iN'*iN);       % estimate R, noise on current
    end
end

% estimate the pdf
BinWidth=2;               % set the parameters for the HIST call
X=[0:BinWidth:1100];
for k=1:3
    R1pdf(:,k)=hist(R1(:,k),X)/BinWidth/NRepeat;   % estimated pdf
    R2pdf(:,k)=hist(R2(:,k),X)/BinWidth/NRepeat;   % estimated pdf
end


% plot the results
FigNum=4;

figure(FigNum)
clf

for k=1:3
    subplot(3,1,k)
    plot(X,R2pdf(:,k),'k',X,R1pdf(:,k),'k')
    axis([900 1100 0 0.04])
    DG_SetFontSize(12)
    DG_SetTraceWidth(2,1,FigNum,k)
    DG_SetTraceWidth(1,2,FigNum,k)
    DG_SetTraceTint(30,1,FigNum,k)
    
end
subplot(3,1,1),title('pdf of the estimates')
subplot(3,1,3),xlabel('R')

% export the figures
% DG_Init4PDF(gcf, 4.5, 4.5/1.618/2);        % fixing the size, half standard height
% DG_MakePDF('FigIdent3Bis.pdf', gcf);

