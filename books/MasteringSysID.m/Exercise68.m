% Chapter 4 Exercise 68
% Identification under feedback conditions using time domain methods
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
N=5000;                 % length of the processed data records
                        % N=5000 or N=50000 
NTrans=500;             % eliminate transients of the simulation
NRep=100;               % number of repeated simulations

fMax=0.4;                % max frequency of the generation filter for random excitation
fMaxPlot=N*min(fMax*1.5,0.5);
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi;   %freq grid

% feedback sysem
OrderG=1;    % order of the system
[bFF,aFF]=cheby1(OrderG,20,0.5);bFF(1)=0;bFF=bFF*10;     % Feed Forward-system
                                                         % unit feedback
br=bFF;ar=aFF+bFF;                                       % r --> y
G0=freqz(bFF,aFF,w);                                     % FRF of the exact feedback system

% G generator
[bGen,aGen]=butter(3,2*fMax);           % filter used to generate the rand reference signal

% noise system
stdNoiseAll=[0.1 0.5 1];                % standard deviation of disturbing noise
OrderNoise=1;                           % BJ noise model
[bNoise,aNoise]=butter(OrderNoise,2*fMax/3,'low');     % v=He noise model
bv=conv(aFF,bNoise);av=conv(aFF+bFF,aNoise);           % v --> y

% Generate and process the data with ARX, OE, and Box-Jenkins 
for r=1:length(stdNoiseAll)             % loop over the noise levels
    stdNoise=stdNoiseAll(r);            % select noise level
    
% system response        
    for s=NRep:-1:1                     % loop over de realizations
    [r s]                           % follow the simulation
    r0=randn(N+NTrans,1);r0=filter(bGen,aGen,r0);  % reference signal
    v=randn(N+NTrans,1)*stdNoise;                  % driving disturbing noise 
    y=filter(br,ar,r0)+filter(bv,av,v);            % response of the system
    u=r0-y;                                        % input of the FF   
    u(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate transients of the simulation

    % three different estimators
    Moe = oe([y u],[OrderG OrderG 1],'maxiter',200,'tol',1e-6,'lim',0);     % outpute-error method
    Marx = arx([y u],[OrderG OrderG 1],'tol',1e-6,'lim',0);                 % arx-method
    Mbj= bj([y u],[OrderG OrderNoise OrderNoise OrderG 1],'maxiter',200,'tol',1e-6,'lim',0);
                                                                            % Box-Jenkins method
    % extract result  
    % Box Jenkins processing
    [b,a]=tfdata(Mbj,'v');                  % import the estimated coefficients          
    G=freqz(b,a,w);                         % calculate estimated FRF
    Gbj(:,s)=G;                             % store the result of simulation j
    ebj(:,s)=abs(G0-G);                     % complex error
    
    % OE processing
    [b,a]=tfdata(Moe,'v');
    G=freqz(b,a,w);
    Goe(:,s)=G;
    eoe(:,s)=G0-G;
       
    % ARX processing
    [b,a]=tfdata(Marx,'v');
    G=freqz(b,a,w);
    Garx(:,s)=G;
    earx(:,s)=G0-G;
       
    end  % end loop over the repetitions
  
    % calculate bias
    arxBias{r}=G0(:)-mean(Garx,2);
    oeBias{r}=G0(:)-mean(Goe,2);
    bjBias{r}=G0(:)-mean(Gbj,2);
    
    %Calculate mean square error
    Ebj{r}=sqrt(mean(abs(ebj).^2,2));
    Eoe{r}=sqrt(mean(abs(eoe).^2,2));
    Earx{r}=sqrt(mean(abs(earx).^2,2));
    
end   % loop over the noise levels

% plot the results the results
FigNum=15;       % plot the RMS-errors
figure(FigNum),clf
counter=0;
f=w/2/pi;

Noise{1}='std 0.1';Noise{2}='std 0.5';Noise{3}='std 1';
for r=1:3   % loop over the noise filters
    counter=counter+1;
    subplot(1,3,counter)
    plot(f,db(G0),'k',f,db(Earx{r}),'k',f,db(Eoe{r}),'k',f,db(Ebj{r}),'k')
    axis([0 fMax -80 20])
    if r==1;ylabel('Amplitude');end
    title(['Noise filter  ' Noise{r}])
    xlabel('Frequency');
    DG_SetFontSize(12,FigNum,r)
       
end

counter=0   % set colors of the figures
for r=1:3   % loop over the noise filters
        counter=counter+1;
        DG_SetTraceWidth(1,'*',FigNum,counter)
        DG_SetTraceTint(60,'*',FigNum,counter)
        DG_SetTraceTint(100,1,FigNum,counter)    % G0
        DG_SetTraceStyle('--',2,FigNum,counter)  % arx model
        DG_SetTraceWidth(1,1,FigNum,counter)
        DG_SetTraceStyle(':',3,FigNum,counter)  % oe model
        DG_SetTraceTint(100,3,FigNum,counter)
        DG_SetTraceWidth(1.25,3,FigNum,counter)
        DG_SetTraceWidth(0.5,4,FigNum,counter)   % bj
        DG_SetTraceTint(100,4,FigNum,counter)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigFeedbackARXandOEandBJRMS.pdf', gcf);   

FigNum=16;       % plot the bias errors
figure(FigNum),clf
counter=0;
f=w/2/pi;

Noise{1}='std 0.1';Noise{2}='std 0.5';Noise{3}='std 1';
for r=1:3   % loop over the noise filters
    counter=counter+1;
    subplot(3,1,counter)
    plot(f,db(G0),'k',f,db(arxBias{r}),'k',f,db(oeBias{r}),'k',f,db(bjBias{r}),'k')
    axis([0 fMax -80 20])
    ylabel('Amplitude (dB)');       
    title(['Noise filter  ' Noise{r}])
    xlabel('Frequency');
    DG_SetFontSize(12,FigNum,r)
end

counter=0   % set colors of the figures
for r=1:3   % loop over the noise filters
        counter=counter+1;
        DG_SetTraceWidth(1,'*',FigNum,counter)
        DG_SetTraceTint(60,'*',FigNum,counter)
        DG_SetTraceTint(100,1,FigNum,counter)    % G0
        DG_SetTraceStyle('--',2,FigNum,counter)  % arx model
        DG_SetTraceWidth(1,1,FigNum,counter)
        DG_SetTraceStyle(':',3,FigNum,counter)  % oe model
        DG_SetTraceTint(100,3,FigNum,counter)
        DG_SetTraceWidth(1.25,3,FigNum,counter)
        DG_SetTraceWidth(0.5,4,FigNum,counter)   % bj
        DG_SetTraceTint(100,4,FigNum,counter)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigFeedbackARXandOEandBJBias5000.pdf', gcf);     
  