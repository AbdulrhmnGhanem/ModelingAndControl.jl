% Chapter 1 Exercise 9
% Identification in the presence of outliers
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010


clear all
NData=100                   % number of measurements in 1 experiment
NRepeat=10000               % Number of repeated experiments
iMax=0.01                   % maximum amplitude of the current
stdU=1                      % standard deviation of the noise
R0=1000                     % exact value of the resistor

Rsearch=R0*[0.8:0.001:1.2]; % search points for a scanned search

for k=NRepeat:-1:1          % repeated experiments
    if rem(k,1000)==0,k,end % display the actual experiment number
    i=iMax*rand(NData,1);   % generate current, uniformly distributed
    
    n=randn(NData,1).^2;    % CHI-Square noise 1 degree of freedom
    nMean=1;nMedian=0.455;  % mean value = 1; std=2; median=0.455 (from Table)
  
    uLap=R0*i+n-nMedian;    % generate Laplace estimator output --> take out median
    uGauss=R0*i+n-nMean;    % generate Gaussian estimator output --> take out mean
    
    % estimate LS
    RLSGauss(k)=uGauss'*i/(i'*i);
    RLSLap(k)=uLap'*i/(i'*i);
    
    % estimate Least Absolute values by simple scan
    for p=1:length(Rsearch);K(p)=sum(abs(uGauss-Rsearch(p)*i));end     % Gaussian data
    [Kmin,P]=min(K);
    RLAVGauss(k)=Rsearch(P);  
    
    for p=1:length(Rsearch);K(p)=sum(abs(uLap-Rsearch(p)*i));end     % Laplace data
    [Kmin,P]=min(K);
    RLAVLap(k)=Rsearch(P); 
    
    
end

% calculate the pdf
N=length(RLSGauss);     % # of points
BinWidth=1;             % set the parameters for the HIST call
X=[800:BinWidth:1200];
pdfLSG=hist(RLSGauss',X)/BinWidth/N;
pdfLAVG=hist(RLAVGauss',X)/BinWidth/N;
pdfLSL=hist(RLSLap',X)/BinWidth/N;
pdfLAVL=hist(RLAVLap',X)/BinWidth/N;


% plot results
FigNum=11;
figure(FigNum)

subplot(1,2,1),plot(X,pdfLSG,'k',X,pdfLAVG,'k')
axis([800 1200 0 0.05]);
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetTraceTint(30,2,FigNum,1)
xlabel('R'),ylabel('pdf')

title('Calibration with mean value')
'mean and standard deviation, LS, calibr. mean',disp(mean(RLSGauss)),disp(std(RLSGauss))
'mean and standard deviation, LAV, calibr. mean',disp(mean(RLAVGauss)),disp(std(RLAVGauss))


subplot(1,2,2)
plot(X,pdfLSL,'k',X,pdfLAVL,'k')
title('Calibration with median')
'mean and standard deviation, LS, calibr. median',disp(mean(RLSLap)),disp(std(RLSLap))
'mean and standard deviation, LS, calibr. median',disp(mean(RLAVLap)),disp(std(RLAVLap))
axis([800 1200 0 0.05]);
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetTraceTint(30,2,FigNum,2)
xlabel('R'),ylabel('pdf')

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent6.pdf', gcf);


