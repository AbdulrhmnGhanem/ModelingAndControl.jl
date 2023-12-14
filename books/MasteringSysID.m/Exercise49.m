% Chapter 3 Exercise 49
% Direct measurement of the FRF under feedback conditions
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 29 November 2010

clear all

% define the system and calculate the transfer functions FRF
[bFF,aFF]=cheby1(2,20,0.5);bFF(1)=0;bFF=bFF*2;     % plant (FF-system), controler C=1
br=bFF;ar=aFF+bFF;                                 % r --> y
bv=aFF;av=aFF+bFF;                                 % v --> y

% define the input variables
NPeriod=1024;               % number of data points per block or period
M=256;                      % number of blocks to be processed
N=NPeriod*M;                % number of data points in the measurement

fs=128                      % sampling frequency 
AmplAll=[0 0.2 0.4 0.8];    % RMS-values of the r-signal
snv=0.2;                    % std. dev. process noise
NTrans=1024                 % eliminate the transient effects simulation

%
%%
%%%
% Part I
%
% FRF measurement with random noise excitation  
% Diff window, 
% M block processed, 
%%%
%%
%

% generate excitation signal
for s=length(AmplAll):-1:1      % loop over the amplitudes of the excitation
    Ampl=AmplAll(s);
    r=randn(N+NTrans,1)*Ampl;   % generate white random noise reference signal
    v=randn(N+NTrans,1)*snv;    % disturbing noise to be added
    
% system response    
    y=filter(br,ar,r)+filter(bv,av,v);  % response of the closed loop system
    u=r-y;                              % input of the FF-branche
    
    y(1:NTrans)=[];u(1:NTrans)=[];      % eliminate transients simulation
       
 % process the data
    u=reshape(u,NPeriod,M);             % one block per collumn
    y=reshape(y,NPeriod,M);
    
    
    Lines=[1:NPeriod/2];                % fft lines used in the analysis
    Y=fft(y)/(NPeriod/2);               % fft analysis
    U=fft(u)/(NPeriod/2);
    Y=diff(Y);U=diff(U);                % apply diff window
    YD=Y(Lines,:);UD=U(Lines,:);
    
    f=(Lines-1)'/NPeriod*fs;            % frequencies for rectangular window
    fD=f+0.5/NPeriod*fs;                % frequencies for diff window
    
    UU=mean(abs(UD(:,1:M)).^2,2);       % SUU with diff window
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2);  % SYU with diff window
    GD(:,s)=YU./UU;                     % FRF estimate
end                                     % end loop over MAll        
   
    GD0=freqz(bFF,aFF,fD/fs*2*pi);      % exact FRF
 
    
%
%%
%%%
% Part II
%
% FRF measurement with multisine excitation as function of SNR
% Rectangular window
% M realizationa
%%%
%%
%



% generate excitation signal

    Lines=[2:NPeriod/2-1];              % fft lines used in the analysis
    f=(Lines-1)'/NPeriod*fs;            % frequencies of the multisine
    
 for s=length(AmplAll):-1:1             % loop over the amplitudes
    Ampl=AmplAll(s);
    r=RandMulti(Lines,NPeriod,N+NTrans)*Ampl;r=r(:);     % generate the random phase multisine
    v=randn(N+NTrans,1)*snv;            % generate the disturbing noise
    
% system response    
    y=filter(br,ar,r)+filter(bv,av,v);  % response of the closed loop system
    u=r-y;                              % input of the FF-branche
    y(1:NTrans)=[];u(1:NTrans)=[];      % eliminate transients
   
 % process the data
    u=reshape(u,NPeriod,M);             % 1 period per collumn
    y=reshape(y,NPeriod,M);
  
    Y=fft(y)/(NPeriod/2);               % fft analysis
    U=fft(u)/(NPeriod/2);
    Y=Y(Lines,:);U=U(Lines,:);  
    
    Um=mean(U,2);                       % mean over the periods 
    Ym=mean(Y,2);
    GMulti(:,s)=Ym./Um;                 % estimated FRF
end      % end loop over MAll
    
GMulti0=freqz(bFF,aFF,f/fs*2*pi);       % exact value FRF at multisine frequencies
    
% plot the results 
FigNum=1;
figure(FigNum)
clf
 
for k=1:4  % amplitudes
    subplot(2,2,k)
    plot(fD,db(GD(:,k)),'k',f,db(GMulti(:,k)),'k',fD,db(GD0),'k')
    ylabel('Amplitude (dB)')
    xlabel('Frequency (Hz)')
    title(['Amplitude ',num2str(AmplAll(k))])
    axis([0 60 -40 20])
    DG_SetFontSize(14)
    DG_SetTraceTint(70,1,FigNum,k)
    DG_SetTraceTint(30,3,FigNum,k)
    DG_SetTraceWidth(2,'*',FigNum,k)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation331.pdf', gcf);

 
 FigNum=2;
 figure(FigNum)
 for k=1:4  % phase
    subplot(2,2,k)
    plot(fD,180/pi*angle(GD(:,k)),'.k',f,180/pi*angle(GMulti(:,k)),'.k',fD,180/pi*angle(GD0),'k')
    ylabel('Phase (deg)')
    if k>2,xlabel('Frequency (Hz)'),end
    title(['Amplitude ',num2str(AmplAll(k))])
    axis([0 60 -180 180])
    DG_SetFontSize(12)
    DG_SetTraceTint(70,1,FigNum,k)
    DG_SetTraceTint(30,3,FigNum,k)
    DG_SetTraceWidth(2,'*',FigNum,k)
 end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation332.pdf', gcf);
  