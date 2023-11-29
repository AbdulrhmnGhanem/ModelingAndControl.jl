% Chapter 1 Exercise 10
% Influence of the number of parameters on the model uncertainty 
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all

Ndata=1000          % number of data points
Nrepeat=100000      % number of repeated experiments

a0=1;b0=0;          % system y=a0 u+b0


u1=linspace(0,1,Ndata);u1=u1(:);    % generate a uniformly distributed input in [0,1]

y10=a0*u1+b0;       % calculate the exact system output


for r=Nrepeat:-1:1              % make Nrepeat experiments
	noise=randn(size(y10));
	y1=y10+noise;               % generate disturbed output
	
	K2=[u1 ones(size(u1))];     % estimate slope and offset
    K1=u1;                      % estimate only slope, offset fixed to zero
	
	Est1=K1\y1;   % This is the numerical stable calculation of (K1'K1)^-1*K1'y1
	Est2=K2\y1;
	
	Est1All(r,:)=Est1';
	Est2All(r,:)=Est2';
end

% calculate pdf
N=length(Est1All(:,1));     % # of points
BinWidth=0.01;              % set the parameters for HIST call
X=[0.5:BinWidth:1.5];
pdfa1=hist(Est1All(:,1),X)/BinWidth/N;
pdfa2=hist(Est2All(:,1),X)/BinWidth/N;

% plot the results
FigNum=12;

figure(FigNum),clf
subplot(1,2,1)
plot(X,pdfa1,'k',X,pdfa2,'k')
title('Distribution of the estimated slope')
xlabel('Slope'),ylabel('pdf')
axis([0.5  1.5 0 10]);
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetTraceTint(30,2,FigNum,1)

subplot(1,2,2)
u1=linspace(0,2,40);u1=u1(:);
u2=u1+(u1(2)-u1(1))/2;
clear ym1All1 ym2All1
for r=100:-1:1
    ym2All(:,r)=(Est2All(r,1)*u1'+Est2All(r,2)-a0*u1')';
    u2Plot(:,r)=u1;
    ym1All(:,r)=(Est1All(r,1)*u2'-a0*u2')';
    u1Plot(:,r)=u2;
end

plot(u1Plot(:),ym1All(:),'.k',u2Plot(:),ym2All(:),'.k')
axis([0 2 -1 1])
xlabel('Time')
ylabel('Error')

DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetTraceTint(30,2,FigNum,2)
xlabel('Time (s)'),title('Error modeled water level')
DG_SetFontSize(12)

% calculate the covariance matrix

C1=cov(Est1All);
C2=cov(Est2All);
disp('mean value slope (1 parameter model'),disp(mean(Est1All))
disp('mena value slope (2 parameter model'),disp(mean(Est2All(:,1)))
disp('Standard deviation slope (2 parameter model)'),disp(sqrt(C2(1,1)))
disp('Standard deviation slope (1 parameter model'),disp(sqrt(C1))


% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIdent7.pdf', gcf);
