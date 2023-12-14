% Chapter 1 Exercise 7
% Chacterizing a 2-dimensional parameter estimate
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
Ndata=100           % number of data points
Nrepeat=1000        % number of repeated experiments

a0=0.1;b0=0;        % true system y0=a0 u + b0      y= y0 + v
stdNoise=1;         % standard deviation of the disturbing noise


u1=linspace(-3,3,Ndata);u1=u1(:);   % first set of input data
u2=linspace(2,5,Ndata);u2=u2(:);    % second set of input data

y10=a0*u1+b0;       % create first set of data
y20=a0*u2+b0;       % create second set of data

for r=Nrepeat:-1:1     % loop over the repeated experiments
	noise=stdNoise*randn(size(y10));    % create disturbing noise
	y1=y10+noise;
	y2=y20+noise;
	
	K1=[u1 ones(size(u1))];             % create K-matrix first data set
	K2=[u2 ones(size(u2))];             % create K-matrix second data set
	
	Est1=K1\y1;   % This is the numerical stable calculation of (K1'K1)^-1*K1'y1
	Est2=K2\y2;
	
	Est1All(r,:)=Est1';
	Est2All(r,:)=Est2';
end


% create the first plot
FigNum=8;

figure(FigNum)
clf
plot(Est2All(:,1),Est2All(:,2),'.k',Est1All(:,1),Est1All(:,2),'.k')
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetTraceTint(30,1,FigNum,1)
xlabel('Slope')
ylabel('Offset')

% Export the first plot
%DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
%DG_MakePDF('FigIdent8a.pdf', gcf);



% evaluate the models in the input range [0,10]
u1=linspace(0,10,40);u1=u1(:);  % create interleaved data set
u2=u1+(u1(2)-u1(1))/2;


for r=100:-1:1  % evaluate all the estimated models
    ym1All1(:,r)=(Est1All(r,1)*u1'+Est1All(r,2))';
    u1Plot(:,r)=u1;
    ym2All1(:,r)=(Est2All(r,1)*u2'+Est2All(r,2))';
    u2Plot(:,r)=u2;
end

% Create the 2nd plot
FigNum=9;
figure(FigNum)
clf

subplot(2,2,1)      % output of model fitted in [-3,3]
plot(u1Plot,ym1All1,'k')
axis([0 10 -5 5])

DG_SetFontSize(12)
DG_SetTraceWidth(1,'*',FigNum,1)
% DG_SetTraceTint(30,2,FigNum,2)
ylabel('Modeled output'),%xlabel('Time (s)')

subplot(2,2,2)       % output of model fitted in [3,5]
plot(u2Plot,ym2All1,'k')
axis([0 10 -5 5])
DG_SetFontSize(12)
DG_SetTraceWidth(1,'*',FigNum,2)
DG_SetTraceTint(30,'*',FigNum,2)

subplot(2,2,3)      % plot the errors on the modelled output
u1=linspace(0,10,30);u1=u1(:);  % create interleaved data set
u2=u1+(u1(2)-u1(1))/2;

clear u1Plot u2Plot
for r=100:-1:1      % calculate the errors
    e1All1(:,r)=(Est1All(r,1)*u1'+Est1All(r,2))'-(a0*u1'+b0)';
    u1Plot(:,r)=u1;
    e2All1(:,r)=(Est2All(r,1)*u2'+Est2All(r,2))'-(a0*u2'+b0)';
    u2Plot(:,r)=u2;
end

plot(u2Plot(:),e2All1(:),'.k',u1Plot(:),e1All1(:),'.k')
axis([0 10 -1 1])
DG_SetFontSize(12)
xlabel('Time')
ylabel('Error y')

DG_SetTraceWidth(1,'*',FigNum,3)
DG_SetTraceTint(30,1,FigNum,3)
xlabel('Time (s)')

subplot(2,2,4)  % plot the std. dev. of the errors
plot(u1,std(e1All1,0,2),'k',u2,std(e2All1,0,2),'k')
axis([0 10 0 1])
DG_SetFontSize(12)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetTraceTint(30,2,FigNum,4)
xlabel('Time (s)')
ylabel('std(y)')

%DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
%DG_MakePDF('FigIdent8b.pdf', gcf);

%print -dpdf FigIdent8



disp('Covariance matrix experiment 1'),disp(cov(Est1All))
disp('Covariance matrix experiment 2'),disp(cov(Est2All))
disp('Correlation matrix experiment 1'),disp(corrcoef(Est1All))
disp('Correlation matrix experiment 1'),disp(corrcoef(Est2All))

