% Chapter 4 Exercise 73
% Comparison of the time and frequency domain identification under feedback
% conditions, using periodic excitations
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
NPer=1024;                      % period length
M=7;                            % number of processed periods
N=M*NPer;                       % total number of processed data
NRep=100;                       % number of repeated simulations
NTrans=NPer;                    % eliminate transient of the simulation
fMax=0.4;                       % max frequency of the generation filter for the excitation
Lines=1+[1:floor(fMax*NPer)];   % freq. lines to be used for Freq. Dom. estimator ELiS
f=(Lines-1)/NPer;               % frequencies to be used for ELiS

% feedback system
OrderG=2;                       % order of the system
[bFF,aFF]=cheby1(OrderG,20,0.5);bFF(1)=0;bFF=bFF*2;     %feed forward-system
                                                        % unit feedback
br=bFF;ar=aFF+bFF;                                      % r --> y

% G generator
[bGen,aGen]=butter(3,2*fMax);   % filter used to generate the excitation of the system

% G0-FRF for plotting
fMaxPlot=N*min(fMax*1.5,0.5);
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi;   %freq grid
G0=freqz(bFF,aFF,w);

% noise models
stdNoiseAll=[0.1 0.5 1];        % standard deviation of disturbing noise
OrderNoise=2;                   % disturbing noise filter
[bNoise,aNoise]=butter(OrderNoise,2*fMax/3,'low');     % v=He noise model
bv=conv(aFF,bNoise);av=conv(aFF+bFF,aNoise);           % v --> y

% Generate and process the data with Box-Jenkins and with Frequency Domain
% identification method
for r=1:length(stdNoiseAll)     % loop over the noise levels
    stdNoise=stdNoiseAll(r);    % select noise level
% system response    
    
    for s=NRep:-1:1             % loop over de realizations
        [r s]                   % follow the simulation
         
     r0=RandMulti(Lines,NPer,NPer*(M+1));                   % flat random phase multisine
        
     v=randn(1,N+NTrans)*stdNoise;                          % driving disturbing noise

     y=filter(br,ar,r0)+filter(bv,av,v);                    % response of the system
     u=r0-y+randn(size(y))*1e-3;                            % input of the FF

     u(1:NTrans)=[];y(1:NTrans)=[];                         % eliminate transient of the simulation
     
     u=reshape(u,NPer,M);U=fft(u)/sqrt(N);U=U(Lines,:);     % fft input
     y=reshape(y,NPer,M);Y=fft(y)/sqrt(N);Y=Y(Lines,:);     % fft output
        
     u=mean(u,2);y=mean(y,2);     % use the averaged values in time domain

    % time domain: Box-Jenkins 
    Mbj= bj([y u],[OrderG OrderNoise OrderNoise OrderG 1],'maxiter',1000,'tol',1e-6,'lim',0);
     
    % frequency domain identification 
    Fdat=fiddata(num2cell(Y,1),num2cell(U,1),f(:));         % store the data in fiddata
    Fdat=varanal(Fdat);                                     % add variance analysis

    nb=OrderG;      % order numerator
    na=OrderG;      % order denominator
    figure(1);      % plot for intermediate results
    MFD=elis(Fdat,'z',nb,na,struct('fs',1));                % frequency domain identification
    
    % Extract the results
    % Box Jenkins processing
    [b,a]=tfdata(Mbj,'v');            
    G=freqz(b,a,w);
    Gbj(:,s)=G;
    ebj(:,s)=G0-G;
       
    % Freq Dom processing
    [domain,bFD,aFD]=imppar(MFD);    % processing the estimates
    G=freqz(bFD,aFD,w);
    GFD(:,s)=G;
    eFD(:,s)=(G0-G);
    
end  % end loop over the repetitions


    
    % calculate bias
    bjBias{r}=G0(:)-mean(Gbj,2);
    ELiSBias{r}=G0(:)-mean(GFD,2);
    
    %Calculate mean square error
    Ebj{r}=sqrt(mean(abs(ebj).^2,2));
    EELiS{r}=sqrt(mean(abs(eFD).^2,2));
    
end   % loop over the noise levels

% Plot the results
FigNum=22;
figure(FigNum)
clf
counter=0;
f=w/2/pi;                       % frequencies to be plotted

Noise{1}='std 0.01';Noise{2}='std 0.1';Noise{3}='std 1';
for r=1:3   % loop over the noise filters
    counter=counter+1;
    subplot(1,3,counter)
    plot(f,db(G0),'k',f,db(Ebj{r}),'k',f,db(EELiS{r}),'k')
    axis([0 fMax -80 20])
    if r==1;ylabel('Amplitude (dB)');end
    title(['Noise filter  ' Noise{r}])
    DG_SetFontSize(12,FigNum,r)
end
counter=0   % set colors of the figures
for r=1:3   % loop over the noise filters
    counter=counter+1;
    DG_SetTraceWidth(2,'*',FigNum,counter)
    DG_SetTraceTint(60,'*',FigNum,counter)
    DG_SetTraceTint(100,1,FigNum,counter)    % G0

    DG_SetTraceWidth(2.5,2,FigNum,counter)   % bj
    DG_SetTraceTint(40,2,FigNum,counter)

    DG_SetTraceWidth(1,3,FigNum,counter)   % ELiS
    DG_SetTraceTint(100,3,FigNum,counter)
    DG_SetTraceStyle(':',3,FigNum,counter)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigFeedbackELiSvsBJRMS.pdf', gcf);  

FigNum=23;
figure(FigNum),clf
counter=0;
f=w/2/pi;

Noise{1}='std 0.1';Noise{2}='std 0.5';Noise{3}='std 1';
for r=1:3   % loop over the noise filters
    counter=counter+1;
    subplot(1,3,counter)
    plot(f,db(G0),'k',f,db(bjBias{r}),'k',f,db(ELiSBias{r}),'k')
    axis([0 fMax -80 20])
    if r==1;ylabel('Amplitude (dB)');end
    title(['Noise filter  ' Noise{r}])
    DG_SetFontSize(18,FigNum,r)
end
counter=0   % set colors of the figures
for r=1:3   % loop over the noise filters
   counter=counter+1;
    DG_SetTraceWidth(2,'*',FigNum,counter)
    DG_SetTraceTint(60,'*',FigNum,counter)
    DG_SetTraceTint(100,1,FigNum,counter)    % G0

    DG_SetTraceWidth(2.5,2,FigNum,counter)   % bj
    DG_SetTraceTint(40,2,FigNum,counter)

    DG_SetTraceWidth(1,3,FigNum,counter)   % ELiS
    DG_SetTraceTint(100,3,FigNum,counter)
    DG_SetTraceStyle(':',3,FigNum,counter)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigFeedbackELiSvsBJBias.pdf', gcf);         

FigNum=24;
figure(FigNum)
clf

plot(f,db(Gbj),'k')
DG_SetFontSize(7,FigNum)
xlabel('Frequency');ylabel('Amplitude (dB)')
DG_SetTraceWidth(1.5,'*',FigNum)
DG_SetTraceTint(30,'*',FigNum)

hold on
plot(f,db(GFD),'k')
hold off       

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigFeedbackELiSvsBJOutliers.pdf', gcf);      