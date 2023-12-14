% Chapter 4 Exercise 75
% Using the time domain identification toolbox from Matlab
%
% The authors thank Lennart Ljung from Linkoping University, Sweden for his 
%     detailed advices to set up this exercise
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 

clear all

% define the input variables
N=1024*7;                               % length of the data set
NTrans=1000;                            % transient points
NRep=10;                                % number of repeated simulations
stdNoise=0.1;                           % standard deviation of disturbing noise
OrderG=2;                               % order of the plant
OrderNoise=1;                           % order of the noise filter

% plant model
[b0,a0]=cheby1(OrderG,10,2*0.25);b0(2)=b0(2)*1.3; % plant to be estimated

% noise model
[b,a]=butter(OrderNoise,2*0.2); b=b+0.001*a;% noise model
bNoise=b;   
aNoise=a; 

% Generator filter
[bGen,aGen]=butter(3,0.3);

% generate the signals
u0=filter(bGen,aGen,randn(1,N+NTrans));u0=u0/std(u0);
y0=filter(b0,a0,u0);                    % exact output
yNoise=filter(bNoise,aNoise,randn(1,N+NTrans));
yNoise=stdNoise*yNoise;                 % disturbing noise
y=y0+yNoise;                            % disturbed output

u0(1:NTrans)=[];y(1:NTrans)=[];y0(1:NTrans)=[];         % eliminate transients

% prepare the data for the time domain toolbox
tData=iddata(y(:),u0(:),1); % create time domain data object
tDataEst=tData(1:N/2);      % estimation set --> to estimate the model
tDataVal=tData(N/2+1:end);  % validation set --> to validate the model
y0Val=y0(N/2+1:end);        % validation data without noise --> to verify the result

% Step 1: take a look at the data
FigNum=1
figure(FigNum)
plot(tDataEst)
title('Estimation data')
shg

% step 2: make an initial guess of the delay. Is a delay present?
nk=delayest(tDataEst)

% Step 3: try a default state space model 
%               - using the initial guess for the delay
%               - using the default order 2
m=pem(tDataEst,'nk',nk);    %  order of the model is selected by the routine

% Step 4: try also order 1, 2, and 3
m1=pem(tDataEst,1,'nk',nk);    % 1st order model
m2=pem(tDataEst,2,'nk',nk);    % 1st order model
m3=pem(tDataEst,3,'nk',nk);    % 3th order model

% Step 5: compare the models on the validation data
FigNum=FigNum+1;
figure(FigNum)
compare(tDataVal,m1,m2,m3)
% m1 is poor, m2 and m3 looks slightly better than m2 
%                (this can vary from one realization to the other)

% From the fit values (see plot) it m3 seems to be the best --> try also m4
m4=pem(tDataEst,4,'nk',nk);    % 3th order model

FigNum=FigNum+1;
figure(FigNum)
compare(tDataVal,m1,m2,m3,m4)
%  m3 turns out to be slightly better 
%                 (this can vary from one realization to the other)
% --> we pick model order 3 to continue

% Study the residuals of model m3

FigNum=FigNum+1;
figure(FigNum)
resid(tDataVal,m3)
% m3 passes 
%    - the whiteness test
%    - the cross-correlation with the inputs

% Compute the spectral analysis of the estimate

FigNum=FigNum+1;
figure(FigNum)
gs=spafdr(tDataEst)
bode(gs,m3,'sd',3,'fill',{0.05,pi,100})

% there is a good agreement between the model and the spectral estimate
% note that you can zoom on the plot

%
% Study of the plant model
%
% study of the pole/zero plot
pzmap(m3)   

% a pole zero cancelation becomes visible
% --> add the uncertainty regions and zoom
FigNum=FigNum+1;
figure(FigNum)
pzmap(m3,'sd',3)
% The uncertainty bounds overlap. This is a strong indication that one
% pole/zero can be eliminated --> model order 2 is also a good candidate model

% 
% study of the noise model
%
FigNum=FigNum+1;
figure(FigNum)
pzmap(m3('n'))   

% a pole zero cancelation becomes visible
% --> add the uncertainty regions and zoom
pzmap(m3('n'),'sd',3)

% The uncertainty bounds overlap for the two complex poles/zeros. 
% This is a strong indication that these pole/zero can be eliminated 
%    --> a first order noise model seems to be enough

% 
% Estimate a Box-Jenkins model:  3th degree plant model, 1st degree noise model
%
FigNum=FigNum+1;
figure(FigNum)
mbj2=bj(tDataVal,[3 1 1 2 0]);
mbj3=bj(tDataVal,[4 1 1 3 0]);
compare(tDataVal,m3,mbj2,mbj3)

% both models result in the same fit --> both models are acceptable
FigNum=FigNum+1
resid(tDataVal,mbj2);     % residue test of bj2 model

% compare the simulated output to y0
ym3=predict(tDataVal,m3,inf);    % calculate the simulated output m3
ybj2=predict(tDataVal,mbj2,inf);
y3=ym3.y;ybj=ybj2.y;             % extract the simulated outputs
sqrt(sum((y3-y0Val(:)).^2/length(y3)))
sqrt(sum((ybj-y0Val(:)).^2/length(y3)))


% prepare plots for book
FigNum=1; %  estimation data
figure(FigNum)    
DG_SetTraceColor('black','*',FigNum,'*')
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('Lennart1.pdf', gcf); 

FigNum=2; % results model m1, m2, m3
figure(FigNum)    
DG_SetTraceColor('black','*',FigNum,'*')
DG_SetTraceTint(100,1,FigNum,1)
DG_SetTraceTint(70,2,FigNum,1)
DG_SetTraceTint(50,3,FigNum,1)
DG_SetTraceTint(30,4,FigNum,1)
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('Lennart2.pdf', gcf); 

FigNum=3; %  % results model m1, m2, m3
figure(FigNum)    
DG_SetTraceColor('black','*',FigNum,'*')
DG_SetTraceTint(100,1,FigNum,1)
DG_SetTraceTint(70,2,FigNum,1)
DG_SetTraceTint(50,3,FigNum,1)
DG_SetTraceTint(30,4,FigNum,1)
DG_SetTraceTint(10,5,FigNum,1)
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('Lennart3.pdf', gcf); 

FigNum=4; %
figure(FigNum)
DG_SetFontSize(14,FigNum,'*')
subplot(2,1,1), title('Correlation function of residuals. Output y1')
xlabel('lag')
subplot(2,1,2), title('Cross corr. function between input u1 and residuals from output y1')
xlabel('lag')
%DG_Init4PDF(gcf, 8);        % fixing the size, half standard height
%DG_MakePDF('Lennart4.pdf', gcf); 

FigNum=6; %
figure(FigNum)   
DG_SetTraceColor('black','*',FigNum,'*')
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('Lennart6.pdf', gcf); 

FigNum=7; % pole zero plot noise model
figure(FigNum)   
DG_SetTraceColor('black','*',FigNum,'*')
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('Lennart7.pdf', gcf); 

