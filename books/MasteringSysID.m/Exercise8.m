% Chapter 1 Exercise 8
% Dependence of the optimal cost function on the distribution of the disturbing noise 
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
NData=100                  % number of measurements in 1 experiment
NRepeat=100000             % Number of repeated experiments
iMax=0.01                  % maximum amplitude of the current
stdU=1                     % standard deviation of the noise
R0=1000                    % exact value of the resistor

Rsearch=R0*[0.9:0.001:1.1];     % search points for a scanned search

for k=NRepeat:-1:1              % loop over the repeated experiments
    if rem(k,1000)==0,disp(k),end  % display the actual experiment number
    i=iMax*rand(NData,1);       % generate current, uniformly distributed
    
    x=rand(NData,1);            % generate Laplace noise
    nLap=zeros(size(x));
    nLap(x<=0.5) = log(2*x(x<=0.5))/sqrt(2)*stdU; 
    nLap(x>0.5) =  - log(2*(1-x(x>0.5)))/sqrt(2)*stdU;
    
    nGauss=randn(size(x))*stdU; % generate Gaussian noise
    
    uLap=R0*i+nLap;             % generate Laplace disturbed output
    uGauss=R0*i+nGauss;         % generate Gaussian disturbed output
    
    % estimate LS
    RLSGauss(k)=uGauss'*i/(i'*i);
    RLSLap(k)=uLap'*i/(i'*i);
    
    % estimate Least Absolute values by simple scan
    for p=1:length(Rsearch);K(p)=sum(abs(uGauss-Rsearch(p)*i));end  % Gaussian data
    [Kmin,P]=min(K);
    RLAVGauss(k)=Rsearch(P);  
    
    for p=1:length(Rsearch);K(p)=sum(abs(uLap-Rsearch(p)*i));end    % Laplace data
    [Kmin,P]=min(K);
    RLAVLap(k)=Rsearch(P); 
    
    
end

% calculate pdf
N=length(RLSGauss);     % # of points
BinWidth=1;             % set paramaters HIST call
X=[900:BinWidth:1100];
pdfLSG=hist(RLSGauss',X)/BinWidth/N;
pdfLAVG=hist(RLAVGauss',X)/BinWidth/N;
pdfLSL=hist(RLSLap',X)/BinWidth/N;
pdfLAVL=hist(RLAVLap',X)/BinWidth/N;

% plot the results
FigNum=10
figure(FigNum)

subplot(1,2,1),
plot(X,pdfLSG,'k',X,pdfLAVG,'k')
axis([900 1100 0 0.05]);
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum)
DG_SetTraceTint(30,2,FigNum)
title('Gaussian noise')
xlabel('R')
'mean and standard deviation Gaussion noise, LS',disp(mean(RLSGauss)),disp(std(RLSGauss))
'mean and standard deviation Gaussion noise, LAV',disp(mean(RLAVGauss)),disp(std(RLAVGauss))
subplot(1,2,2)
plot(X,pdfLSL,'k',X,pdfLAVL,'k')
axis([900 1100 0 0.05]);
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum)
DG_SetTraceTint(30,2,FigNum)
title('Laplace noise')
xlabel('R')
'mean and standard deviation Laplace noise, LS',disp(mean(RLSLap)),disp(std(RLSLap))
'mean and standard deviation Laplace noise, LAV',disp(mean(RLAVLap)),disp(std(RLAVLap))

% Export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent5.pdf', gcf);

