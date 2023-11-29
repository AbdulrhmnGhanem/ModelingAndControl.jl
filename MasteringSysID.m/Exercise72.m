% Chapter 4 Exercise 72
% Emphasizing a frequency band using nonparametric noise models and
% periodic excitations
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
NPer=1024;                          %  period length
M=7;                                % number of processed periods
MTrans=1;MAll=M+MTrans;             % number of periods for transient elimination of the simulation
NTrans=NPer*MTrans;                 % transient points
N=M*NPer;                           % total number of processed data
NRep=100;                          % number of repeated simulations
stdNoise=0.1;                       % standard deviation of disturbing noise
OrderG=2;                           % order of the plant
OrderNoise=2;                       % BJ noise model

fGen=0.25;                          % bandwidth of the generator signal (fs=1)
OrderPF=3;                          % order of the prefilter 
fmax=0.1;                           % frequency band of the prefilter

fMaxPlot=N*0.5;                     % frequency range used for plotting
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi; % freq grid

[bGen,aGen]=butter(3,2*fGen);bGen(2)=0.9*bGen(2);   % generator filter
Ggen=abs(freqz(bGen,aGen,w));                       % FRF generator filter

[b0,a0]=cheby1(OrderG,5,2*0.08);b0(2)=b0(2)*1.3;    % plant G0 to be estimated
G0=freqz(b0,a0,w);                                  % FRF of the plant

[b,a]=butter(OrderNoise,2*0.2); b=b+0.1*a;          % disturbing noise filter
bNoise=b;   
aNoise=a;
Gnoise=stdNoise*abs(freqz(bNoise,aNoise,w));        % FRF of the disturbing noise filter

[bf,af]=butter(OrderPF,fmax*2);bf=bf+0.01*af;       % create lowpss prefilter
GPF=freqz(bf,af,w);                                 % FRF prefilter

LinesCut=floor([2:1:fmax*NPer]);                    % FFT lines used with prefiltering
fCut=(LinesCut-1)/NPer;                             % frequencies to be used for ELiS
Lines=1+[1:NPer/2-1];                               % FFT lines used if no prefiltering
f=(Lines-1)/NPer;                                   % frequencies to be used for ELiS

SNR=db(Ggen./Gnoise);   % signal-to-noise-ratio

for s=NRep:-1:1      % loop over de realizations
    [s]
      
    % generate the signals, no prefiltering
    u0=RandMulti(Lines,NPer,NPer*MAll);             % flat random phase multisine
    u0=filter(bGen,aGen,u0);                        % filtered input signal
    y0=filter(b0,a0,u0);                            % exact output
    yNoise=filter(bNoise,aNoise,randn(1,N+NTrans));yNoise=stdNoise*yNoise;
                                                    % disturbing noise
    y=y0+yNoise;                                    % disturbed ouptut
    
    % generate the prefiltered signals
    u0PF=filter(bf,af,u0);yPF=filter(bf,af,y);      % prefilter the input and output
    
    % process the signals
    u0(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate transients
    
    % no prefiltering
    u0=reshape(u0,NPer,M);y=reshape(y,NPer,M);      % 1 period per collumn
    U0=fft(u0)/sqrt(N);U0Full=U0(Lines,:);          % FFT input
    Y=fft(y)/sqrt(N);YFull=Y(Lines,:);              % FFT output

    Fdat=fiddata(num2cell(YFull,1),num2cell(U0Full,1),f(:));    
                                                    % store the data in fiddata
    Fdat=varanal(Fdat);                             % add variance analysis
   
    % Cut frequency band
    YCut=Y(LinesCut,:);U0Cut=U0(LinesCut,:);
    FdatCut=fiddata(num2cell(YCut,1),num2cell(U0Cut,1),fCut(:));   
                                                    % store the data in fiddata
    FdatCut=varanal(FdatCut);                       % add variance analysis
 
    % prefiltered signals
    u0PF(1:NTrans)=[];yPF(1:NTrans)=[];             % eliminate transients simulation
    u0PF=reshape(u0PF,NPer,M);yPF=reshape(yPF,NPer,M);
                                                    % 1 period per collumn
    U0PF=fft(u0PF)/sqrt(N);U0PF=U0PF(Lines,:);      % FFT input
    YPF=fft(yPF)/sqrt(N);YPF=YPF(Lines,:);          % FFT output

    FdatPF=fiddata(num2cell(YPF,1),num2cell(U0PF,1),f(:));    
                                                    % store the data in fiddata
    FdatPF=varanal(FdatPF);                         % add variance analysis

    % frequency domain identification 
    nb=OrderG;      % order numerator
    na=OrderG;      % order denominator
    figure(1)       % plot for intermediate results
    MFD=elis(Fdat,'z',nb,na,struct('fs',1));        % no prefiltering
    MFDCut=elis(FdatCut,'z',nb,na,struct('fs',1));  % cut the frequency band
    MFDPF=elis(FdatPF,'z',nb,na,struct('fs',1));    % prefiltered data

    % Extract the results
    [domain,bFD,aFD]=imppar(MFD);                   % extract the estimates
    G=freqz(bFD,aFD,w);                             % FRF
    GFD(:,s)=G;
    eFD(:,s)=(G0-G);
    
    [domain,bFDPF,aFDPF]=imppar(MFDPF);    
    G=freqz(bFDPF,aFDPF,w);
    GFDPF(:,s)=G;
    eFDPF(:,s)=(G0-G);
    
    [domain,bFDCut,aFDCut]=imppar(MFDCut);   
    G=freqz(bFDCut,aFDCut,w);
    GFDCut(:,s)=G;
    eFDCut(:,s)=(G0-G);
    
end  % end loop over the repetitions
    
%Calculate mean square error
EFD=sqrt(mean(abs(eFD).^2,2));
EFDCut=sqrt(mean(abs(eFDCut).^2,2));
EFDPF=sqrt(mean(abs(eFDPF).^2,2));

% Figures for book
FigNum=21;
figure(FigNum)
clf

f1=w/2/pi;  % plot frequencies
plot(f1,db(G0),'k',f1,db(GPF),'k',f1,db(EFD),'k',f1,db(EFDCut),'k',f1,db(EFDPF),'k')
axis([0 0.5 -80 10])
ylabel('Amplitude');
xlabel('Frequency');
DG_SetFontSize(7,FigNum)
DG_SetTraceTint(100,1,FigNum)     
DG_SetTraceWidth(1.25,1,FigNum)    
DG_SetTraceTint(100,2,FigNum)     
DG_SetTraceWidth(0.5,2,FigNum)       
DG_SetTraceTint(100,3,FigNum)     
DG_SetTraceWidth(0.5,3,FigNum)      
DG_SetTraceTint(50,4,FigNum)     
DG_SetTraceWidth(1,4,FigNum)      
DG_SetTraceTint(100,5,FigNum)     
DG_SetTraceWidth(0.5,5,FigNum)      
DG_SetTraceStyle(':',5,FigNum)


% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigELiSPrefiltering.pdf', gcf); 
      