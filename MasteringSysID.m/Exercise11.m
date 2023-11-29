% Chapter 1 Exercise 11
% Model selection using the AIC-criterion 
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
NEst=1000                       % length estimation record
NValid=10000                    % length validation record

stdyAll=[0.5  0.05];            % std. dev. distrubing noise at the output
nbMax=100                       % maximum order of the FIR-model
Order=[1:nbMax];                % orders to be scanned

[b,a]=cheby1(3,0.5,2*[0.15 0.3]);   % the coefficients of the true IIR system


% generate the test data with random input
u=randn(1,NEst);                % exact input
y0=filter(b,a,u);               % exact output
yNoise=randn(size(y0));         % noise disturbance to be added to the output


% generate the validation data
uValid=randn(1,NValid);             % exact input
yValid0=filter(b,a,uValid);         % exact output
yValidNoise=randn(size(yValid0));   % disturbed output

% process the data

for qSTD=length(stdyAll):-1:1   % run over the different noise levels
    
    y=y0+yNoise*stdyAll(qSTD);                  % add disturbing noise to test data
    yValid=yValid0+yValidNoise*stdyAll(qSTD);   % add disturbing noise to validation data

	% generate the maximum size K matrix
	KAll=zeros(NEst,nbMax);
	for s=1:NEst
            p=[1:min(s,nbMax+1)];
            q=s+1-p;
            KAll(s,p)=u(q);
	end
	
	% estimate the different orders
	
	for r=length(Order):-1:1
        nb=Order(r);
        K=KAll(:,1:nb);       
        bEst=K\y'; % make least squares estimate
	
		% check the cost on the estimation data
		yEst=filter(bEst,1,u);
		eEst=y-yEst;
		
		
		Est(r,qSTD)=sum(eEst.^2)/NEst/stdyAll(qSTD)^2;  % normalized estimation cost
		AICRel(r,qSTD)=Est(r,qSTD)*(1+2*nb/NEst);       % AIC-criterion
		
		
		% check the cost on the validation data
		yModel=filter(bEst,1,uValid);
		eValid=yValid-yModel;
		Valid(r,qSTD)=sum(eValid.^2)/NValid/stdyAll(qSTD)^2;    % normalized validation cost
        
        % noise less validation data
       	eNoiseLess=yValid0-yModel;
		VNL(r,qSTD)=sqrt(sum(eNoiseLess.^2)/NValid/stdyAll(qSTD)^2); % normalized validation cost	
	end


end

FigNum=1;
figure(FigNum),clf


    
subplot(2,2,1)
plot(Order,Est(:,1),'-k',Order,Valid(:,1),'-k',Order, AICRel(:,1),'-k')
axis([1 100 0.7 1.2])
DG_SetFontSize(12)
ylabel('Cost')
title('Noisy data')
DG_SetTraceTint(30,3,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)

subplot(2,2,2)
plot(Order,VNL(:,1),'-k')
axis([1 100 0 1])
DG_SetFontSize(12)
ylabel('Normalized rms')
title('Noiseless data')
DG_SetTraceWidth(2,'*',FigNum,2)
    
    
subplot(2,2,3)
plot(Order,Est(:,2),'-k',Order,Valid(:,2),'-k',Order, AICRel(:,2),'-k')
axis([1 100 0.7 1.2])
DG_SetFontSize(12)
xlabel('Order')
ylabel('Cost')
DG_SetTraceTint(30,3,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)

subplot(2,2,4)
plot(Order,VNL(:,2),'-k')
axis([1 100 0 1])
DG_SetFontSize(12)
xlabel('Order')
ylabel('Normalized rms')
DG_SetTraceWidth(2,'*',FigNum,4)

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigAIC.pdf', gcf);
    
