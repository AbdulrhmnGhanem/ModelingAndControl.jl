% Chapter 3 Exercise 45
% Impact of the block (period) length on the uncertainty
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

NPeriodAll=[128  512 2048 8192];    % number of data points per block (period)
M=64;                               % number of blocks to be processed

fs=128                              % sampling frequency
Ampl=0.33                           % RMS-value of the multi sine
sny=0.02;                           % std. dev. noise on the output
NTrans=1024                         % eliminate transient effects of the simulation

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
for r=length(NPeriodAll):-1:1       % loop over the block length
    
    NPeriod=NPeriodAll(r);
    N=NPeriod*M;                    % Number of data points in the measurement
   
    u=randn(N+NTrans,1)*Ampl;
    
% system response    
    y0=filter(b,a,u);                % response of the system
    yNoise=sny*randn(size(y0));      % noise to be added to the output
    y=y0+yNoise;                     % noisy output
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
    YD=Y(Lines,:);UD=U(Lines,:);Y0D=Y0(Lines,:);    % select the frequencies to be processed
    
    f=(Lines-1)'/NPeriod*fs;          % frequencies for rectangular window
    fD{r}=f+0.5/NPeriod*fs;           % frequencies for diff window
    
    
    UU=mean(abs(UD(:,1:M)).^2,2);     % SUU process with diff window
    YY=mean(abs(YD(:,1:M)).^2,2);     % SYY process with diff window
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2);      % SYU
    Y0U=mean(Y0D(:,1:M).*conj(UD(:,1:M)),2);    % SY0U
    GD{r}=YU./UU;                     % frf estimate for noisy data
    G0D=Y0U./UU;                      % frf estimate for noiseless data
                              
    Coh=abs(YU).^2 ./(abs(UU).*abs(YY));  % estimate the theoretical variance from the last realization
    GD0Var{r}=abs(GD{r}).^2.*(1-Coh)./Coh/M;GD0Std{r}=sqrt(GD0Var{r});   % Theoretic std. dev.
    
    GD0{r}=freqz(b,a,fD{r}/fs*2*pi);   % exact  FRF
       
end   % end loop over NPeriodAll
    

%
%%
%%%
% FRF measurement with multisine excitation.  Variance analysis for a fixed M
% Rectangular window
% noise added to the output
%%%
%%
%


% generate multisine excitation signal    
 clear f   
    
 for r=length(NPeriodAll):-1:1          % loop over the period length
     
    NPeriod=NPeriodAll(r); 
    Lines=[2:NPeriod/2-1];              % fft-lines to be plotted
    N=NPeriod*M;                        % total record length
   
    u=RandMulti(Lines,NPeriod,N+NTrans)*Ampl;     % generate the multisine excitation

 % system response    

    y0=filter(b,a,u);                   % response of the system
    yNoise=sny*randn(size(y0));         % noise to be added to the output
    y=y0+yNoise;                        % add the noise to the output
    y(1:NTrans)=[];u(1:NTrans)=[];      % eliminate transients
   
 % process the data

    u=reshape(u,NPeriod,M);             % one period per collumn
    y=reshape(y,NPeriod,M);
     
    Y=fft(y)/(NPeriod/2);               % fft analysis
    U=fft(u)/(NPeriod/2);
    Y=Y(Lines,:);U=U(Lines,:);          % select the frequencies
    
    f{r}=(Lines-1)'/NPeriod*fs;         % frequency grid corresponding to the record length NPeriod

    Um=mean(U,2);                       % mean values 
    Ym=mean(Y,2);
    GMulti{r}=Ym./Um;                   % FRF estimate

    GMulti0{r}=freqz(b,a,f{r}/fs*2*pi); % exact FRF
    
    YStd=std(Y,0,2);                    % estimate the theoretical variance of FRF from last realization
    GMulti0Std{r}=abs(GMulti{r}).*abs(YStd./Ym)/sqrt(M);
    
 end                             % end loop over NPeriodAll   

% plot figures for noise and multisine excitation as function NPeriodAll 
FigNum=5
figure(FigNum)   
clf
 
for r=1:4 
subplot(2,2,r)
    plot(fD{r},db(GD{r}),'k',fD{r},db(GD0Std{r}),'k',f{r},db(GMulti0Std{r}),'k')
    ylabel('Amplitude (dB)')
    xlabel('Frequency (Hz)')
    axis([0 60 -60 1])
    title([num2str(NPeriodAll(r)),' Block length'])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,r)
    DG_SetTraceWidth(2,1,FigNum,r)
    DG_SetTraceWidth(2,2,FigNum,r)
end

% export plots    
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation315.pdf', gcf);

 
 
 
 


