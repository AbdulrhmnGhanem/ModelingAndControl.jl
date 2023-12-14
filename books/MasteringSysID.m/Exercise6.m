% Chapter 1 Exercise 6
% Least squares estimation of models that are linear in the parameters
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
N=100;                          % number of data points
nMax=20;                        % max degree of the polynomial model
u0=linspace(0,1,N);u0=u0(:);    % input data
y0=tan(u0*0.9*pi/2);            % nonlinear system
y=y0;                           % no noise added to the output

% setup of the maximum sized K-matrix  
for k=nMax:-1:1
    KAll(:,k)=u0.^(k-1);
end

% run over the model of decreasing degree
for k=nMax:-1:1
    K=KAll(:,[1:k]);
    CondK(k,2)=cond(K);
    CondK(k,1)=cond(K'*K);
    
    % numerical unstable solution
    bEst1=(inv(K'*K))*(K'*y);bEst1=bEst1(:)';
    yMod1=polyval(fliplr(bEst1),u0);
    e1=y0-yMod1;
    e1All(k)=sqrt(sum(e1.^2)/N);    % mean square error
    
    % numerical stable solution
    bEst2=K\y0;bEst2=bEst2(:)';
    yMod2=polyval(fliplr(bEst2),u0);
    e2=y0-yMod2;
    e2All(k)=sqrt(sum(e2.^2)/N);    % mean square error
end


% plot the results
OrderAll=[1:nMax]-1;

FigNum=8;
figure(FigNum),clf

subplot(1,2,1)
semilogy(OrderAll,CondK(:,1),'k',OrderAll,CondK(:,2),'k');
DG_SetFontSize(12)
xlabel('Order')
ylabel('Cond. nr.')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)

subplot(1,2,2)
semilogy(OrderAll,abs(e1All),'k',OrderAll,abs(e2All),'k');
DG_SetFontSize(12)
xlabel('Order')
ylabel('rms error')
DG_SetTraceTint(30,2,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2) 

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent10.pdf', gcf);
 





