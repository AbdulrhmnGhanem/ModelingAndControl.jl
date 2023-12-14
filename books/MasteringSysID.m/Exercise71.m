% Chapter 4 Exercise 71
% Identification in the frequency domain using nonparametric noise models
% and periodic data
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
NPer=1024;                                  % period length
M=7;                                        % number of processed periods
N=M*NPer;                                   % total number of processed data
NTrans=NPer;                                % to eliminate transients of the simulation
NRep=20;                                    % number of repeated simulations
SignalSelect=2;                             % 1 --> periodic noise   2 --> random phase multisine
OrderG=2;                                   % order of the plant
fMax=0.4;                                   % bandwidth of the generator signal (fs=1)
stdNoise=0.1;                               % standard deviation of disturbing noise

Lines=1+[1:floor(fMax*NPer)];               % freq. lines in FFT
f=(Lines-1)/NPer;                           % frequencies to be used for ELiS

fMaxPlot=N*min(fMax*1.5,0.5);               % used for plotting
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi;         %freq grid

[bGen,aGen]=butter(3,2*fMax);               % filter used to generate the excitation of the system
[b0,a0]=cheby1(OrderG,10,2*fMax*0.9);       % plant to be estimated
G0=freqz(b0,a0,w);                          % FRF plant

% define the three noise models  ARX, OE, BJ
bNoise{1}=1;        % ARX noise model
aNoise{1}=a0;
OrderNoise(1)=OrderG;                       % model order for BJ-estimation

bNoise{2}=1;        % OE noise model --> white noise
aNoise{2}=1;
OrderNoise(2)=0;                            % model order for BJ-estimation

OrderNoise(3)=2;    % BJ noise model
[b,a]=butter(OrderNoise(3),2*fMax/3,'low'); % noise model
bNoise{3}=b;   
aNoise{3}=a; 

% generate and process the data in the time domain and in the frequency
% domain
for r=1:3             % loop over the different noise filters (== time domain method ARX OE Box-Jenkins)
     Ggen=abs(freqz(bGen,aGen,w));Gnoise=stdNoise*abs(freqz(bNoise{r},aNoise{r},w)); % generator and noise filter
     SNR{r}=db(Ggen./Gnoise);   % signal-to-noise-ratio
     
    for s=NRep:-1:1     % loop over de realizations
        [r,s]
        u0=RandMulti(Lines,NPer,NPer*(M+1));                    % flat random phase multisine
        u0=filter(bGen,aGen,u0);                                % filtered input
        y0=filter(b0,a0,u0);                                    % exact output
        yNoise=filter(bNoise{r},aNoise{r},randn(1,N+NTrans));yNoise=stdNoise*yNoise;
                                                                % disturbing noise
        y=y0+yNoise;                                            % disturbed output
        u0(1:NTrans)=[];y(1:NTrans)=[];                         % eliminate transients of the simulation
        u0=u0(:);y=y(:);                                        % make colum vectors
        
        u0=reshape(u0,NPer,M);                                  % one block per collumn
        y=reshape(y,NPer,M);
        
        U0=fft(u0)/sqrt(N);U0=U0(Lines,:);                      % FFT input
        Y=fft(y)/sqrt(N);Y=Y(Lines,:);                          % FFT output
        
        u0=mean(u0,2);y=mean(y,2);     % use the averaged values in time domain
        
        % three different estimators in time domain
        switch r   % select the time domain method
            case 1  % ARX
            MTD = arx([y u0],[OrderG OrderG+1 0],'tol',1e-6,'lim',0);
            case 2   % Output error
            MTD = oe([y u0],[OrderG+1 OrderG 0],'maxiter',1000,'tol',1e-6,'lim',0);
            case 3   % Box-Jenkins
            MTD= bj([y u0],[OrderG+1 OrderNoise(r) OrderNoise(r) OrderG 0],'maxiter',1000,'tol',1e-6,'lim',0);
        end
        
        % frequency domain identification
   
        Fdat=fiddata(num2cell(Y,1),num2cell(U0,1),f(:));    % store the data in fiddata
        Fdat=varanal(Fdat);                                 % add variance analysis 
        nb=OrderG;                                          % order numerator
        na=OrderG;                                          % order denominator
        
        figure(1)                                           % figure to show intermediate iteration results
        MFD=elis(Fdat,'z',nb,na,struct('fs',1));
        
        
% extract results
    
    % Time Domain processing
    [bTD,aTD]=tfdata(MTD,'v');            
    G=freqz(bTD,aTD,w);
    GTD(:,s)=G;
    eTD(:,s)=(G0-G);

    
    % Freq Dom processing
    [domain,bFD,aFD]=imppar(MFD);    % processing the estimates
    G=freqz(bFD,aFD,w);
    GFD(:,s)=G;
    eFD(:,s)=(G0-G);
 
    end  % end loop over the repetitions
    
    %Calculate mean square error
    ETD{r}=sqrt(mean(abs(eTD).^2,2));
    EFD{r}=sqrt(mean(abs(eFD).^2,2));
end      % end loop over the different noise models

% plot the results
FigNum=20;
figure(FigNum),clf
counter=0;
f=w/2/pi;                               % frequencies for plots

Noise{1}='1/A (ARX)';Noise{2}='white (OE)';Noise{3}='C/D (BJ)';
for r=1:3   % loop over the noise filters
    counter=counter+1;
    subplot(1,3,counter)
    plot(f,db(G0),'k',f,db(ETD{r}),'k',f,db(EFD{r}),'k',f,-SNR{r},'k')
    axis([0 fMax -80 1])
    ylabel='Amplitude';
    xlabel='Frequency';
    title(['Noise filter  ' Noise{r}])  
    DG_SetFontSize(12,FigNum,r)
end
counter=0   
for r=1:3   % loop over the noise filter       
    counter=counter+1;
    DG_SetTraceWidth(2,'*',FigNum,counter)
    DG_SetTraceTint(100,1,FigNum,counter)     
    DG_SetTraceWidth(2.5,2,FigNum,counter)      
    DG_SetTraceTint(40,2,FigNum,counter)
    DG_SetTraceStyle(':',3,FigNum,counter)   
    DG_SetTraceTint(100,3,FigNum,counter)
    DG_SetTraceWidth(2.5,3,FigNum,counter)
    DG_SetTraceWidth(1,4,FigNum,counter)    
    DG_SetTraceTint(50,4,FigNum,counter)        
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigELiSvsARXandBJandOE.pdf', gcf);  
