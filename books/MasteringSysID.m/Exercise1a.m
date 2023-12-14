% Chapter 1 Exercise 1.a
% Least squares estimation of the value of a resistor
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 23 November 2010


clear all


% define the input variables 

NRepeat=100                 % Number of repeated experiments
NData=[10 100 1000 10000]   % number of measurements in 1 experiment
iMax=0.01                   % maximum amplitude of the current
stdU=1                      % std. dev. disturbing noise
R0=1000                     % exact value resistor

for r=1:length(NData)       % run of the different data lengths
	for k=1:NRepeat         % run over the repeated experiments
        i=iMax*2*(rand(NData(r),1)-0.5);    % generate random current
        u=R0*i+stdU*randn(size(i));   % generate disturbed output
        R(k,r)=u'*i/(i'*i); % calculate the LS-estimate
	end
end

% plot the results

FigNum=1;

figure(FigNum)   % plot the repeated estimates for different values of N'
clf
K=[1:NRepeat]';


subplot(2,2,1)
plot(K,ones(size(K))*R0,'k',K,R(:,1),'.k')
axis([0 NRepeat 0.9*R0 1.1*R0])
xlabel('Exp. Nr.'),ylabel('R(10)')
DG_SetTraceTint(30,1,FigNum,1)
DG_SetTraceWidth(3,1,FigNum,1)

subplot(2,2,2)
plot(K,ones(size(K))*R0,'k',K,R(:,2),'.k')
axis([0 NRepeat 0.9*R0 1.1*R0])
xlabel('Exp. Nr.'),ylabel('R(100)')
DG_SetTraceTint(30,1,FigNum,2)
DG_SetTraceWidth(3,1,FigNum,2)

subplot(2,2,3)
plot(K,ones(size(K))*R0,'k',K,R(:,3),'.k')
axis([0 NRepeat 0.9*R0 1.1*R0])
xlabel('Exp. Nr.'),ylabel('R(1000)')
DG_SetTraceTint(30,1,FigNum,3)
DG_SetTraceWidth(3,1,FigNum,3)

subplot(2,2,4)
plot(K,ones(size(K))*R0,'k',K,R(:,4),'.k')
axis([0 NRepeat 0.9*R0 1.1*R0])
xlabel('Exp. Nr.'),ylabel('R(10 000)')
DG_SetTraceTint(30,1,FigNum,4)
DG_SetTraceWidth(3,1,FigNum,4)


DG_SetFontSize(14,FigNum,'*')

% export the figures
%  DG_Init4PDF(gcf, 4.5);                                   % fixing the size
%  DG_MakePDF('FigIdent1.pdf', gcf);

