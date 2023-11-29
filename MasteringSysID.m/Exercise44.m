% Chapter 3 Exercise 44
% Analysis of the noise errors in FRF measurements
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables
NPeriod=1024;  % number of data points per block
M=16;          % number of blocks to be processed
NRepeat=10;    % number of repeated experiments to measure the sample variance
fs=128         % sampling frequency s.t. 1 period = 1 s
Ampl=0.33      % RMS-value of the excitation
sny=0.02;      % std. dev. noise= sny*2*median(y0)
NTrans=1024    % to eliminate transient effects


[b,a]=cheby1(2,10,0.2);           % test system

N=NPeriod*M;            % Number of data points in the measurement

%
%%
%%%
% Part I
%
% FRF measurement with random noise excitation.  Variance analysis for a fixed M
% diff window
% noise added to the output
%%%
%%
%

% generate excitation signal
for r=NRepeat:-1:1  
   
    u=randn(N+NTrans,1)*Ampl;
    
% system response    

    y0=filter(b,a,u);                   % response of the system
    yNoise=sny*randn(size(y0));         % noise to be added to the output
    y=y0+yNoise;                        % add disturbing noise to the output
    y0(1:NTrans)=[];y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients
     
 % process the data
 
    u=reshape(u,NPeriod,M);             % 1 block per collumn
    y=reshape(y,NPeriod,M);
    y0=reshape(y0,NPeriod,M);
    
    Lines=[1:NPeriod/2];                % frequencies to be processed
          
    Y=fft(y)/(NPeriod/2);               % fft-analysis
    Y0=fft(y0)/(NPeriod/2);
    U=fft(u)/(NPeriod/2);
    Y=diff(Y);U=diff(U);Y0=diff(Y0);    % apply diff window
    YD=Y(Lines,:);UD=U(Lines,:);Y0D=Y0(Lines,:);
    
    f=(Lines-1)'/NPeriod*fs;            % frequencies fft
    fD=f+0.5/NPeriod*fs;                % frequencies diff-window analysis
    
    
    UU=mean(abs(UD(:,1:M)).^2,2);           % SUU
    YY=mean(abs(YD(:,1:M)).^2,2);           % SYY
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2);  % SYU
    Y0U=mean(Y0D(:,1:M).*conj(UD(:,1:M)),2);% SY0U
    GD(:,r)=YU./UU;G0D=Y0U./UU;
end                                     % end loop over MAll

    Coh=abs(YU).^2 ./(abs(UU).*abs(YY));  % estimate the theoretical variance from the last realization
    GD0Var=abs(GD(:,end)).^2.*(1-Coh)./Coh/M;GD0Std=sqrt(GD0Var);   % Theoretic std. dev.
    
    GD0=freqz(b,a,fD/fs*2*pi);              % G0 at the diff frequencies
    GDM=mean(GD,2);                         % mean value over all realizations
    GDStd=std(GD,0,2);                      % std. dev.
    

%
%%
%%%
% Part II
%
% FRF measurement with multisine excitation.  Variance analysis for a fixed M
% Rectangular window
% noise added to the output
%%%
%%
%


% generate excitation signal

    Lines=[2:NPeriod/2-1];
    
    
 for r=NRepeat:-1:1
   
    u=RandMulti(Lines,NPeriod,N+NTrans)*Ampl;     % generate the excitation

 % system response    

    y0=filter(b,a,u);               % response of the system
    yNoise=sny*randn(size(y0));     % noise to be added to the output
    y=y0+yNoise;                    % add noise to the output
    y(1:NTrans)=[];u(1:NTrans)=[];  % eliminate transients
   
 % process the data

    u=reshape(u,NPeriod,M);         % 1 period per collumn
    y=reshape(y,NPeriod,M);
     
    Y=fft(y)/(NPeriod/2);           % spectral analysis
    U=fft(u)/(NPeriod/2);
    Y=Y(Lines,:);U=U(Lines,:);
    
    f=(Lines-1)'/NPeriod*fs;         % freq. vector


    Um=mean(U,2);                    % mean value 
    Ym=mean(Y,2);
    GMulti(:,r)=Ym./Um;              % FRF estimate
end       % end loop over MAll

    GMulti0=freqz(b,a,f/fs*2*pi);    % exact value
    GMultiM=mean(GMulti,2);          % mean value over all realizations
    GMultiStd=std(GMulti,0,2);       % std. dev.
    
    YStd=std(Y,0,2);                 % variance estiomation of the output noise 
    GMulti0Std=abs(YStd./Um)/sqrt(M);% estimate the theoretical variance of FRF from last realization
    
% plot figures noise and multisine  different M 
FigNum=4;
figure(FigNum)
clf
 
subplot(2,2,1)
    plot(fD,db(GDM-GD0),'k.',fD,db(GDStd/sqrt(NRepeat)),'k',fD,db(GD0),'k')
    ylabel('Amplitude (dB)')
    xlabel('f (Hz)')
    title('Noise excitation')
     axis([0 60 -60 1])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,1,FigNum,1)
    DG_SetTraceWidth(2,'*',FigNum,1)
    
subplot(2,2,2)
    plot(f,db(GMultiM-GMulti0),'k.',f,db(GMultiStd/sqrt(NRepeat)),'k',f,db(GMulti0),'k')
    ylabel('Amplitude (dB)')
    xlabel('f (Hz)')
     axis([0 60 -60 1])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,1,FigNum,2)
    DG_SetTraceWidth(2,'*',FigNum,2)
    title('Multisine excitation')
    
subplot(2,2,3)
    plot(fD,db(GD0Std),'k',fD,db(GDStd),'k',fD,db(GD0),'k')
    ylabel('Amplitude (dB)')
    xlabel('f (Hz)')
      axis([0 60 -60 1])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,1,FigNum,3)
     DG_SetTraceWidth(2,1,FigNum,3)
    DG_SetTraceWidth(2,3,FigNum,3)
   
    
subplot(2,2,4)
    plot(f,db(GMulti0Std),'k',f,db(GMultiStd),'k',f,db(GMulti0),'k')
    ylabel('Amplitude (dB)')
    xlabel('f (Hz)')
      axis([0 60 -60 1])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,1,FigNum,4)
    DG_SetTraceWidth(2,1,FigNum,4)
    DG_SetTraceWidth(2,3,FigNum,4)
 
% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation314.pdf', gcf);

 
 
 
