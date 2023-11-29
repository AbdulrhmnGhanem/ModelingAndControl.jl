% Chapter 3 Exercise 47
% Measuring the FRF in the presence of input and output noise: analysis of
% the errors
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
NRepeat=50;    % number of repeated experiments to measure the sample variance

fs=128         % sampling frequency 
Ampl=1         % RMS-value of the excitation
sny=0.02;      % std. dev. noise output
snu=0.5        % std. dev. noise input
NTrans=1024    % eliminate transient effects simulation

[b,a]=cheby1(2,10,0.2);           % test system


%
%%
%%%
% Part I
%
% FRF measurement with random noise excitation.  Variance analysis for a fixed M
% diff window with varying block length
% noise added to the output
%%%
%%
%

% generate excitation signal
for r=NRepeat:-1:1                   % loop over the realizations
    
    N=NPeriod*M;            % Number of data points in the measurement
    u0=randn(N+NTrans,1)*Ampl;
    
% system response    
    y0=filter(b,a,u0);                % response of the system
    yNoise=sny*randn(size(y0));       % noise to be added to the output
    uNoise=snu*randn(size(y0));       % noise to be added to the output
    
    u=u0+uNoise;                      % noisy input
    y=y0+yNoise;                      % noisy output   
    y0(1:NTrans)=[];y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients
       
 % process the data
    u=reshape(u,NPeriod,M);           % one block per collumn
    y=reshape(y,NPeriod,M);
    y0=reshape(y0,NPeriod,M);
    
    Lines=[1:NPeriod/2];              % fft lines used in the analysis
                                
    Y=fft(y)/(NPeriod/2);             % fft analysis
    Y0=fft(y0)/(NPeriod/2);
    U=fft(u)/(NPeriod/2);
    Y=diff(Y);U=diff(U);Y0=diff(Y0);  % apply diff window
    YD=Y(Lines,:);UD=U(Lines,:);Y0D=Y0(Lines,:);  % select the frequencies to be processed
    
    f=(Lines-1)'/NPeriod*fs;          % frequencies for rectangular window
    fD=f+0.5/NPeriod*fs;              % frequencies for diff window
    
    
    UU=mean(abs(UD(:,1:M)).^2,2);    % SUU process with diff window
    YY=mean(abs(YD(:,1:M)).^2,2);    % SUY process with diff window
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2);      % SYU
    Y0U=mean(Y0D(:,1:M).*conj(UD(:,1:M)),2);    % SUU
    GD(:,r)=YU./UU;G0D=Y0U./UU;       % frf estimate for noisy data
end                                   % end loop over MAll

    Coh=abs(YU).^2 ./(abs(UU).*abs(YY));  % estimate the theoretical variance from the last realization
    GD0Var=abs(GD(:,end)).^2.*(1-Coh)./Coh/M/NRepeat;GD0Std=sqrt(GD0Var);   % Theoretic std. dev.
    
    GD0=freqz(b,a,fD/fs*2*pi);        % exact FRF for diff window frequencies
    GDM=mean(GD,2);                   % mean value over all realizations
    GDStd=std(GD,0,2)/sqrt(NRepeat);  % std. dev.
    

%
%%
%%%
% FRF measurement with multisine excitation.  Variance analysis for a fixed M
% Rectangular window
% noise added to the input and output
%%%
%%
%


% generate excitation signal

 Lines=[2:NPeriod/2-1];
    
    
 for r=NRepeat:-1:1                  % loop over the realizations
   
     u0=RandMulti(Lines,NPeriod,N+NTrans)*Ampl;     % generate the excitation

 % system response    
    y0=filter(b,a,u0);                % response of the system
    yNoise=sny*randn(size(y0));       % noise to be added to the output
    uNoise=snu*randn(size(y0));       % noise to be added to the output
    
    u=u0+uNoise;                      % noisy input
    y=y0+yNoise;                      % noisy output
    
    y(1:NTrans)=[];u(1:NTrans)=[];    % eliminate transients
   
 % process the data
    u=reshape(u,NPeriod,M);           % on period per collumn
    y=reshape(y,NPeriod,M);
                                    
    Y=fft(y)/(NPeriod/2);             % fft analysis
    U=fft(u)/(NPeriod/2);
    Y=Y(Lines,:);U=U(Lines,:);        % slect the frequencies
    
    f=(Lines-1)'/NPeriod*fs;          % frequencies

    Um=mean(U,2);                     % mean over the periods 
    Ym=mean(Y,2);
    GMulti(:,r)=Ym./Um;               % FRF estimate
end      % end loop over NRepeat

GMulti0=freqz(b,a,f/fs*2*pi);         % exact FRF at multisine frequencies      
GMultiM=mean(GMulti,2);               % mean value over the realizations
GMultiStd=std(GMulti,0,2)/sqrt(NRepeat);  % std. dev.

YVar=std(Y,0,2).^2;            % estimate the theoretical variance of FRF from last realization
UVar=std(U,0,2).^2;
clear sYU

% calculate the variance for multisine
for k=length(Lines):-1:1;
    sYU(k,1)=(Y(k,:)-Ym(k))*(U(k,:)-Um(k))'/(M-1);   % covariance
end

GMulti0Var=abs(GMulti(:,end)).^2 .* (YVar./abs(Ym).^2+UVar./abs(Um).^2-2*real(sYU./(Ym.*conj(Um))))/M/NRepeat;
GMulti0Std=sqrt(GMulti0Var);
    

% plot the figures noise and multisine  M=1 and M=16

FigNum=4
figure(FigNum)
clf
 
subplot(2,2,1)
    plot(fD,db(GDM-GD0),'k.',fD,db(GDStd),'k',fD,db(GD0),'k')
    ylabel('Amplitude (dB)')
    xlabel('Frequency (Hz)')
    axis([0 60 -80 1])
     DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,1)
    DG_SetTraceWidth(2,'*',FigNum,1)
    title('Noise excitation')
    
subplot(2,2,2)
    plot(f,db(GMultiM-GMulti0),'k.',f,db(GMultiStd),'k',f,db(GMulti0),'k')
      xlabel('Frequency (Hz)')
    axis([0 60 -80 1])
     DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,2)
    DG_SetTraceWidth(2,'*',FigNum,2)
    title('Multisine excitation')
    
subplot(2,2,3)
    plot(fD,db(GD0Std),'.k',fD,db(GDStd),'k',fD,db(GD0),'k')
      xlabel('Frequency (Hz)')
      ylabel('Amplitude (dB)')
    axis([0 60 -80 1])
     DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,3)
    DG_SetTraceWidth(2,'*',FigNum,3)
    
subplot(2,2,4)
    plot(f,db(GMulti0Std),'.k',f,db(GMultiStd),'k',f,db(GMulti0),'k')
      xlabel('Frequency (Hz)')
      ylabel('Amplitude (dB)')
    axis([0 60 -80 1])
     DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,4)
    DG_SetTraceWidth(2,'*',FigNum,4)

% export the figure
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation324.pdf', gcf);



 
