% Chapter 4 Exercise 67
% Identification in the time domain using parametric noise models
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
N=5000;                                     % length of the processed data records
NTrans=500;                                 % eliminate transients in simulation
NRep=20;                                    % number of repeated simulations

OrderG=2;                                   % order of the plant
fMax=0.4;                                   % bandwidth of the generator signal (fs=1)
[bGen,aGen]=butter(3,2*fMax);               % filter used to generate the excitation of the system
[b0,a0]=cheby1(OrderG,10,2*fMax*0.9);       % plant to be estimated
stdNoise=0.1;                               % standard deviation of disturbing noise

% G0
fMaxPlot=N*min(fMax*1.5,0.5);
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi;         %freq grid
G0=freqz(b0,a0,w);

% noise models
bNoise{1}=1;                                % ARX noise model
aNoise{1}=a0;
OrderNoise(1)=OrderG;                       % model order for BJ noise model

bNoise{2}=1;                                % OE noise model --> white noise
aNoise{2}=1;
OrderNoise(2)=0;

OrderNoise(3)=2;                            % BJ noise model
[b,a]=butter(OrderNoise(3),2*fMax/3,'low'); % noise model
bNoise{3}=b;   
aNoise{3}=a; 


% Simulations: generation and estimation of the data
for r=1:3                                   % loop over the noise filters
     Ggen=abs(freqz(bGen,aGen,w));Gnoise=stdNoise*abs(freqz(bNoise{r},aNoise{r},w)); % generator and noise filter
     SNR{r}=db(Ggen./Gnoise);               % signal-to-noise-ratio
    for s=NRep:-1:1                         % loop over de realizations
        [r,s]                               % follow the simulation
        u0=randn(N+NTrans,1);u0=filter(bGen,aGen,u0);   % filtered random input
        y0=filter(b0,a0,u0);                % exact output
        yNoise=filter(bNoise{r},aNoise{r},randn(N+NTrans,1));
        yNoise=stdNoise*yNoise;             % scaled disturbing noise
        y=y0+yNoise;                        % add disturbing noise to exact output
        u0(1:NTrans)=[];y(1:NTrans)=[];     % eliminate transients of the simulation
        
        % three different estimators
        Moe = oe([y u0],[OrderG+1 OrderG 0],'maxiter',1000,'tol',1e-6,'lim',0);     % output error method
        Marx = arx([y u0],[OrderG OrderG+1 0],'tol',1e-6,'lim',0);                  % ARX-mehtod
        Mbj= bj([y u0],[OrderG+1 OrderNoise(r) OrderNoise(r) OrderG 0],'maxiter',1000,'tol',1e-6,'lim',0);
                                                                                    % Box-Jenkins method
 
% extract results
    
    % Box Jenkins processing
    [b,a]=tfdata(Mbj,'v');                  % import the estimated coefficients          
    G=freqz(b,a,w);                         % calculate estimated FRF
    Gbj(:,s)=G;                             % store the result of simulation j
    ebj(:,s)=abs(G0-G);                     % complex error

    % OE processing
    [b,a]=tfdata(Moe,'v');
    G=freqz(b,a,w);
    Goe(:,s)=G;
    eoe(:,s)=abs(G0-G);
 
    
    % ARX processing
    [b,a]=tfdata(Marx,'v');
    G=freqz(b,a,w);
    Garx(:,s)=G;
    earx(:,s)=abs(G0-G);
    
    end  % end loop over the repetitions
    
    %Calculate mean square error
    Ebj{r}=sqrt(mean(ebj.^2,2));
    Eoe{r}=sqrt(mean(eoe.^2,2));
    Earx{r}=sqrt(mean(earx.^2,2));
end

% plot the results
FigNum=13;
figure(FigNum)
clf

counter=0;
f=w/2/pi;                                    % frequencies for plot

Noise{1}='1/A (ARX)';Noise{2}='white (OE)';Noise{3}='C/D (BJ)';
for r=1:3   % loop over the noise filters
        counter=counter+1;
        subplot(1,3,counter)
        plot(f,db(G0),'k',f,db(Earx{r}),'k',f,db(Eoe{r}),'k',f,db(Ebj{r}),'k',f,-SNR{r},'k')
        axis([0 fMax -80 0])
        if r==1;ylabel('Amplitude (dB)');end
        title(['Noise filter  ' Noise{r}])
        DG_SetFontSize(12,FigNum,r)
end
xlabel('Frequency');
counter=0   % set colors of the figures
for r=1:3   % loop over the noise filters
        counter=counter+1;
        DG_SetTraceWidth(2,'*',FigNum,counter)
        DG_SetTraceTint(60,'*',FigNum,counter)
        DG_SetTraceTint(100,1,FigNum,counter)    % G0
        DG_SetTraceStyle('--',2,FigNum,counter)  % arx model
        DG_SetTraceWidth(2,1,FigNum,counter)
        DG_SetTraceStyle(':',3,FigNum,counter)  % oe model
        DG_SetTraceTint(100,3,FigNum,counter)
        DG_SetTraceWidth(2.5,3,FigNum,counter)
        DG_SetTraceWidth(1,4,FigNum,counter)   % bj
        DG_SetTraceTint(100,4,FigNum,counter)
        DG_SetTraceTint(30,5,FigNum,counter)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigARXandOEandBJ.pdf', gcf);  
   