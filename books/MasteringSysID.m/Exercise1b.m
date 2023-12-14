% Chapter 1 Exercise 1.b
% Analysis of the standard deviation
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 23 November 2010

clear all

% define the input variables 

NRepeat=1000                % Number of repeated experiments
NData=[10 100 1000 10000]   % number of measurements in 1 experiment
iMax=0.01                   % maximum amplitude of the current
stdU=1                      % std. dev. disturbing noise
R0=1000                     % exact value resistor

for r=1:length(NData)       % run over the different data lengths
	for k=1:NRepeat         % run over the repeated experiments
        i=iMax*ones(NData(r),1);    % generate cosntant current
        u=R0*i+stdU*randn(size(i)); % generate disturbed output
        R(k,r)=u'*i/(i'*i); % calculate the LS-estimate
	end
end

stdR=std(R,0,1);     % calculate standard deviation

% Plot the results

FigNum=2;

figure(FigNum)
clf
     
K=logspace(1,log10(NData(end)),100);
loglog(K,stdU./sqrt((K.*iMax^2)),'.k',NData,stdR,'ok') 
title('Standard deviation R(N)')
xlabel('N')
ylabel('std Rest')
DG_SetFontSize(7)
DG_SetTraceTint(30,1,FigNum,1)
DG_SetTraceWidth(3,1,FigNum,1)

% Export the figure
% DG_Init4PDF(gcf, 6);                                   % fixing the size
% DG_MakePDF('FigIdent2.pdf', gcf);

