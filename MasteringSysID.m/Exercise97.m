% Chapter 7 Exercise 97
% Parametric estimation of Gbla
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 8 December 2010


clear all
N=5000;         % length of the processed data records
NTrans=500;     % length transient part to eliminate transients in simulation
NRep=100;       % number of repeated simulations
IterMax=100;    % maximum number of itereations in the Box-Jenkins method

OrderG1=2 ;     % order of the plant G1
OrderG2=2;      % Order of the plant G2
NLpol=[1 1 1];  % coeff. of the static nonlinear polynomial function f
                % x, x^2, x^3
stdNoise=0.05;  % std. dev. of the disturbing white noise added to the output

fMax1=0.25;     % to set the bandwidth of the first systems (fs=1)
fMax2=0.4;      % to set the bandwidth of the first systems (fs=1)
fMaxGen=0.45;   % sets the bandwidth of the generator signal

fPlot=[0:0.001:0.5];            % frequencies used to make the plots
Nnorm=floor(length(fPlot)*0.8); % used to normalize the gain of the estimated transfer function

% define the systems
[bGen,aGen]=butter(3,2*fMaxGen);bGen(end)=bGen(end);    % filter used to generate the excitation of the system
[b10,a10]=butter(OrderG1,2*fMax1);                      % plant to be estimated
[b20,a20]=cheby1(OrderG1,10,2*fMax1);                   % plant to be estimated


% Transfer functions for the plots
wPlot=fPlot*2*pi;           %freq grid
G10=freqz(b10,a10,wPlot);
G20=freqz(b20,a20,wPlot);
GBLA0=G10.*G20;             % GBLA
GGen=freqz(bGen,aGen,wPlot);

% noise models: white Gaussian noise

counter=0;   % count the number of converged runs

for s=NRep:-1:1      % loop over the realizations
      
        % create the data
        u0=randn(N+NTrans,1);u0=filter(bGen,aGen,u0);   % filtered input noise
        p0=filter(b10,a10,u0);                          % first linear system G1
        q0=NLpol(1)*p0+NLpol(2)*p0.^2+NLpol(3)*p0.^3;   % static nonlinearity f
        y0=filter(b20,a20,q0);                          % 2nd linear system G1
        
        y=y0+ stdNoise*randn(size(y0));                 % add white Gaussian noise                                       % no disturbing noise is added
        
        u0(1:NTrans)=[];y(1:NTrans)=[];p0(1:NTrans)=[];q0(1:NTrans)=[];   % eliminate transients
        
%        
% Parametric noise model for BJ 
%
    OrderG=OrderG1+OrderG2;     % order for the BJ-model
    OrderNoise=OrderG+1;        % order that is (arbitrary) selected for the noise model of the nonlinear distortions
                                % it results in an almost white prediction error
    Mbj= bj([y u0],[OrderG+1 OrderNoise OrderNoise OrderG 0],'maxiter',IterMax,'tol',1e-6,'lim',0);
 
% 
 
    if max([Mbj.EstimationInfo.Iterations ])<IterMax   %  check convergence?
                           counter=counter+1;            
        [s counter] 
    
        % Box Jenkins processing
        [b,a]=tfdata(Mbj,'v');  
    
        % extract BLA-model
        G=freqz(b,a,wPlot);
        Gbj(:,counter)=G;
        ebj(:,counter)=abs(GBLA0(:)-Gbj(:,counter));
    
        % extract theoretical std. dev. on amplitude
        [mag,phse,wbj,sdmag,sdphase]=bode(Mbj,wPlot);sdmagAll(:,counter)=sdmag(:);
    
end  % end if loop
    
%        
% Non Parametric noise model with FD-identification
%

% processing with the local polynomial method
U0=fft(u0(:))/sqrt(N);Y=fft(y(:))/sqrt(N);U0(N/2+1:end)=[];Y(N/2+1:end)=[];
Data.U = U0(:).';
Data.Y = Y(:).';            % store spectra as rows
Data.Freq=[0:N/2]'/N;       % store corresponding frequencies

Method.moment=6;            % To have enough degrees of freedom in the result
[CY, Ym, TY, GLocPol1, CvecG2OE,M] = LocalPolyAnalv3(Data,Method); 
YLocPolVar=squeeze(CY.m_nt);% take out the variance
YLocPol=squeeze(Ym.m_nt);

fs=1;fLocPol=[0:N/2-1]'/N*fs;   % create freq. vector Loc. Pol. Method

% Frequency Domain identification
Fdat=fiddata(YLocPol(:),U0,fLocPol,YLocPolVar,0);   % create Frequency Domain Data object
nb=OrderG;                      % order numerator
na=OrderG;                      % order denominator
MFD=elis(Fdat,'z',nb,na,struct('fs',1));
[domain,bFD,aFD]=imppar(MFD);   % processing the estimates
G=freqz(bFD,aFD,wPlot);
GFD(:,s)=G;%/NLGain;
eFD(:,s)=GBLA0(:)-GFD(:,s);     % error on GFD

    
end  % end loop over the repetitions
  
%Calculate mean and std. dev.
stdFD=std(GFD,0,2);    
stdbj=std(Gbj,0,2);
GBJm=mean(Gbj,2);
GFDm=mean(GFD,2);


FigNum=1;                           % compare estimated and theoretical Gbla
figure(FigNum),clf
  
plot(fPlot,db(GBLA0),'k',fPlot,db(GBJm),'k',fPlot,db(GFDm(:)./GBLA0(:)),'--k')
axis([0 0.5 -80 10])
ylabel('Amplitude (dB)')
xlabel('Frequency');
DG_SetFontSize(18,FigNum)        
DG_SetTraceTint(100,1,FigNum)       % GBLA0
DG_SetTraceWidth(2,1,FigNum)         
DG_SetTraceTint(100,2,FigNum)       % GBJm
DG_SetTraceTint(50,3,FigNum)        % db(GFDm(:)./GBLA0(:))
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('BLA_BJ_FD.pdf', gcf);   



FigNum=2                            % study of the variance of parametric Gbla
figure(FigNum),clf
sdmagBJExp=std(abs(Gbj),0,2);       % experimental std. dev. amplitude
sdmagFDExp=std(abs(GFD),0,2);       % experimental std. dev. amplitude
sdmagTh=mean(sdmagAll.^2,2).^0.5;   % averaged mean theoretical value BJ

plot(fPlot,db(GBJm),'k',fPlot,db(sdmagBJExp),'k',fPlot,db(sdmagFDExp),'k',fPlot,db(sdmagTh),'k',fPlot,db(GFDm-GBJm),'k')
axis([0 0.5 -80 10])
ylabel('Amplitude (dB)')
xlabel('Frequency');
DG_SetFontSize(18,FigNum)        
DG_SetTraceTint(100,1,FigNum)       % GBJm
DG_SetTraceWidth(2,1,FigNum)         
DG_SetTraceTint(50,2,FigNum)        % sdmagBJExp
DG_SetTraceTint(50,3,FigNum)        % sdmagFDExp
DG_SetTraceStyle('--',3,FigNum)
DG_SetTraceTint(100,4,FigNum)       % sdmagTh
DG_SetTraceStyle('--',5,FigNum)     % GFDM-GBJ

FigNum=3;                           % analysis of the residuals of last run
figure(FigNum),clf
e=resid(Mbj,[y u0]);
 
FigNum=4;                           % fft-analysis of the prediction errors
figure(FigNum),clf
E=fft(e);
EAmp=db(E(1:N/2)/sqrt(N));fAmp=[0:N/2-1]/N;
Y=fft(y)/sqrt(N);YAmp=db(Y(1:N/2));
plot(fAmp,YAmp,'k',fAmp,EAmp,'k')
axis([0 0.5 -60 20])
ylabel('Amplitude (dB)')
xlabel('Frequency');
DG_SetFontSize(18,FigNum)  
DG_SetTraceTint(50,1,FigNum)    

%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('BLA_BJ_FD_std.pdf', gcf);   

