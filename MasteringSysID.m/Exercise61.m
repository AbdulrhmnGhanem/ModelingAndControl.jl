% Chapter 4 Exercise 61
% Simulation and one-step-ahead prediction
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

%
% This exercise starts from the estimated parameters in Exercise 1
% Run this exercise first
%

disp('run first Exercise 1 of Chapter 4')

% simulation and prediction of the output
ySim=filter(bEst,aEst,u0);                      % simulate the output
eSim=y-ySim;                                    % simulation error
aEstPred=aEst;aEstPred(1)=0;                    % prepare coeff. for prediction filtering
yPred=filter(bEst,1,u0)+filter(-aEstPred,1,y);  % predict the output
ePred=y-yPred;                                  % prediction error
t=[0:N-1];                                      % time vector for plot

% Plot the results
FigNum=3;           % simulation and prediction error
figure(FigNum)
clf
plot(t,[y eSim ePred],'k')
axis([0 50 -1.5 1.5])
DG_SetFontSize(7)
xlabel('time')
DG_SetTraceTint(20,1,FigNum,1)
DG_SetTraceTint(50,3,FigNum,1)
DG_SetTraceTint(100,2,FigNum,1)
DG_SetTraceWidth(1.5,'*',FigNum,1)

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigArxTDSimPredErr.pdf', gcf); 
