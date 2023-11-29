% Chapter 1 Exercise 2
% Study of the asymptotic distribution of an estimate
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 23 November 2010
clear all

% define the input variables

NRepeat=100000    % Number of repeated experiments
NData=[1 2 4 8]   % number of measurements in 1 experiment
iMax=0.01         % maximum amplitude of the current
stdU=1/5          % std. dev. of the disturbing noise
R0=1000           % exact value of the resistor

for r=1:length(NData)               % run over the different data lengths
	for k=NRepeat:-1:1              % run over the repeated experiments
        i=iMax*ones(NData(r),1);    % generate constant current
        uU=R0*i+stdU*(rand(size(i))-0.5)*2*sqrt(3);     % generate uniformly disturbed output
        uN=R0*i+stdU*randn(size(i));                    % generate normaly disturbed output
        RU(k,r)=uU'*i/(i'*i);       % Uniform solution
        RN(k,r)=uN'*i/(i'*i);       % Gaussian solution
	end
end

% calculate the pdf

BinWidth=2;
X=[900:BinWidth:1100];
for k=1:4
    RNpdf(:,k)=hist(RN(:,k),X)/BinWidth/NRepeat;   % estimated pdf
    RUpdf(:,k)=hist(RU(:,k),X)/BinWidth/NRepeat;   % estimated pdf
end


% plot the results
FigNum=3;

figure(FigNum)
clf
for k=1:4
    subplot(4,1,k)
    plot(X,RNpdf(:,k),'k',X,RUpdf(:,k),'k')
    axis([940 1060 0 0.06])
  
    DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,k)
    DG_SetTraceWidth(2,'*',FigNum,k)
end
subplot(4,1,1),title('pdf of the estimates')
subplot(4,1,4),xlabel('R')

% Export the figures
% DG_Init4PDF(gcf, 4.5, 4.5/1.618/2);        % fixing the size, half standard height
% DG_MakePDF('FigIdent3.pdf', gcf);
% print -dpdf FigIdent3


