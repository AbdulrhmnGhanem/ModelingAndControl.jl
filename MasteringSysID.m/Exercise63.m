%  Identification of linear dynamic systems 
%          using an ARX model and no disturbing noise
%
% Chapter 4 Exercise 63
% Sensitivity of the simulation and prediction error to model errors
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010


%
% This exercise starts from the estimated parameters in Exercise 5, Chapter 4
% Run this exercise first
%

disp('run first Exercise 5 of Chapter 4')

% simulation and prediction of the output
ySim=filter(bEst,aEst,u0);                          % calculate the simulated output
eSim=y-ySim;                                        % simulation error
aEstPred=aEst;a(1)=0;                               % prepare coeff. for prediction filtering
yPred=filter(bEst,1,u0)+filter(-aEstPred,1,y)+y;    % calculation of the prediction error
ePred=y-yPred;                                      % prediction error

t=[0:N-1];                                          % time axis for the plot
     
FigNum=2;   % simulation and prediction error
figure(FigNum)
clf

plot(t,[y eSim ePred],'k')
axis([0 100 -1 1])
DG_SetFontSize(7)
xlabel('time')
DG_SetTraceTint(20,1,FigNum,1)   % output
DG_SetTraceTint(50,3,FigNum,1)   % prediction error
DG_SetTraceTint(100,2,FigNum,1)  % simulation error
DG_SetTraceWidth(1.5,'*',FigNum,1)

% export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxModErrPredSim.pdf', gcf);  
% print -dpdf FigArxModErrPredSim

