% Chapter 4 Exercise 66
% One-step-ahead prediction of a noise sequence
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
N=1000;                                         % length data record to be processed
NTrans=1000;                                    % eliminate transient of the simulation
BWAll=[0.01 0.05 0.1 0.2 0.3 0.4];              % noise bandwidth to be checked
t=[1:N];                                        % time vector for plots

% set the order na=5 for first plot;  na=50 for the 2nd plot
na=5;                                           % model order for the arx model (no input)


% create and process the data; plot the results
FigNum=11
figure(FigNum)
clf

FigCounter=0;                                   % initialization of plot number

for k=1:length(BWAll)                           % loop over the bandwidth
BW=BWAll(k);                                    % select the actual bandwidth
[bGen,aGen]=butter(5,2*BW);                     % noise coloring filter
y0=filter(bGen,aGen,randn(1,N+NTrans));         % generate the noise sequency
y0=y0/std(y0);                                  % normalized noise sequence to be modelled
y0ToPred=y0(NTrans+1:NTrans+N);                  % noise to be predicted after elimination transients

% prediction by keeping the value of the last available sample
yEst=y0;
yEst=[0 yEst];                                  % keep the last value
yEst=yEst(NTrans+1:NTrans+N);                   % predicted noise sequence for y0ToPred
e=y0ToPred-yEst;                                 % error on ToPredation set

% prediction using an arx model
yRHS=-y0(na+2:end);                             % right handside of the equations
K=[toeplitz(y0(na+1:end-1),y0(na+1:-1:2))];     % K-matrix, no input available
Theta=K\yRHS(:);                                % least squares solution K Theta = yHRS
aEst=[0 Theta(1:na)'];                          % model parameters for the prediction 
yARX=filter(-aEst,1,y0);                        % arx prediction
eARX=y0-yARX;eARX=eARX(NTrans+1:NTrans+N);      % error predicted value using ARX on y0ToPred
yARX=yARX(NTrans+1:NTrans+N);


 
FigCounter=FigCounter+1;
subplot(3,4,FigCounter)
plot(t,yEst,'k',t,yARX,'k',t,y0ToPred,'k:')
axis([0 floor(N/BW*0.01) -4 4])
title(['Noise  band width ',num2str(BW)])
DG_SetTraceTint(50,1,FigNum,FigCounter)
DG_SetTraceTint(20,2,FigNum,FigCounter)
DG_SetTraceWidth(3,'*',FigNum,FigCounter)
 
FigCounter=FigCounter+1;
subplot(3,4,FigCounter)
plot(t,e,'k',t,eARX,'k')
axis([0 floor(N/BW*0.01) -4 4])
title(['Error band width ',num2str(BW)])
DG_SetTraceTint(50,1,FigNum,FigCounter)
DG_SetTraceTint(20,2,FigNum,FigCounter)
DG_SetTraceWidth(3,'*',FigNum,FigCounter)


RelErr(k,1)=std(e)/std(y0ToPred);
RelErr(k,2)=std(eARX)/std(y0ToPred);

end

DG_SetFontSize(12,FigNum)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigNoisePredna5.pdf', gcf); 


