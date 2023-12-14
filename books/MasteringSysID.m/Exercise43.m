% Chapter 3 Exercise 43
% Measurement of the FRF using a random noise sequence and a random phase
% multisine in the presence of output noise
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

NPeriod=1024;           % number of data points per block
MAll=[1 4 16 64];       % number of blocks to be processed
fs=128                  % sampling frequency s.t. 1 period = 1 s
Ampl=0.33               % RMS-value of the multi sine
sny=0.02;               % std. dev. noise on the output
NTrans=1024             % eliminate transient effects of the simulations

[b,a]=cheby1(2,10,0.2);         % test system

%
% Part 1: white noise excitation and diff window
%

% generate excitation signal
for r=length(MAll):-1:1         % loop over the number of blocks
    M=MAll(r);
    N=NPeriod*M;                % Number of data points in the measurement
    u=randn(N+NTrans,1)*Ampl;   % white noise excitation
    
% system response    
    y0=filter(b,a,u);                 % response of the system
    yNoise=sny*randn(size(y0));       % noise to be added to the output
    y=y0+yNoise;                      % noisy output
    y0(1:NTrans)=[];y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients
      
 % process the data
    u=reshape(u,NPeriod,M);     % one block per collumn
    y=reshape(y,NPeriod,M);
    y0=reshape(y0,NPeriod,M);
    
    Lines=[1:NPeriod/2];        % fft-lines to be used in the processing
    Y=fft(y)/(NPeriod/2);       % fft
    Y0=fft(y0)/(NPeriod/2);
    U=fft(u)/(NPeriod/2);
    Y=diff(Y);U=diff(U);Y0=diff(Y0);                % apply diff window
    YD=Y(Lines,:);UD=U(Lines,:);Y0D=Y0(Lines,:);    % select the lines
    
    f=(Lines-1)'/NPeriod*fs;     % frequency grid 
    fD=f+0.5/NPeriod*fs;         % frequency grid for diff window
    
    
    UU=mean(abs(UD(:,1:M)).^2,2);                   % SUU for diff window
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2);          % SYU for diff window
    Y0U=mean(Y0D(:,1:M).*conj(UD(:,1:M)),2);        % SY0U for diff window
    GD(:,r)=YU./UU;G0D=Y0U./UU;
end                               % end loop over MAll        
   
    GD0=freqz(b,a,fD/fs*2*pi);    % exact FRF for diff window
    
 % figure of a time domain measurement
 FigNum=2;
    figure(FigNum)
    clf
    
    t=[0:NPeriod-1]/fs;   % time axis for plot
    
    clf
    subplot(1,2,1)
    plot(t,u(:,1),'k')
    xlabel('Time (s)')
    ylabel('Input')
    axis([0 8 -1.1 1.1])
    DG_SetFontSize(12)
    DG_SetTraceWidth(2,'*',FigNum,1)

    subplot(1,2,2)
    plot(t,y0(:,1),'k',t,yNoise(1:NPeriod),'k')
    axis([0 8 -1.1 1.1])
    xlabel('Time (s)')
    ylabel('Output')
    DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,2)
    DG_SetTraceWidth(2,'*',FigNum,2)
    
 % DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
 % DG_MakePDF('FigExcitation312.pdf', gcf); 
   
 
    

%
%%
%%%
% part 2: multisine excitation
%
% FRF measurement with multisine excitation.  Show the error for increasing M
% rectangular window
% MAll block processed
% noise added to the output
%%%
%%
%

% generate excitation signal
Lines=[2:NPeriod/2-1];
       
 for r=length(MAll):-1:1                % loop over number of averages
    M=MAll(r); 
    N=NPeriod*M;                        % length record
   
    u=RandMulti(Lines,NPeriod,N+NTrans)*Ampl;     % generate the random phase multisine
 
    % system response    
    y0=filter(b,a,u);                   % response of the system
    yNoise=sny*randn(size(y0));         % nwhite oise to be added to the output
    y=y0+yNoise;                        % noisy output
    y(1:NTrans)=[];u(1:NTrans)=[];      % eliminate transients
   
    % process the data
    u=reshape(u,NPeriod,M);             % 1 period per collumn
    y=reshape(y,NPeriod,M);                          
    Y=fft(y)/(NPeriod/2);               % spectral analysis
    U=fft(u)/(NPeriod/2);
    Y=Y(Lines,:);U=U(Lines,:);          % select the frequencies
    
    f=(Lines-1)'/NPeriod*fs;            % frequency vector

    Um=mean(U,2);                       % mean value over the periods 
    Ym=mean(Y,2);
    GMulti(:,r)=Ym./Um; 
end                             % end loop over MAll
    
GMulti0=freqz(b,a,f/fs*2*pi);           % exact FRF for reference
    
    
% figures noise and multisine  M=1,4,16,64 
FigNum=3;
figure(FigNum);
clf
 
for k=1:4
    subplot(2,2,k)
    plot(fD,db(GD(:,k)-GD0),'k.',f,db(GMulti(:,k)-GMulti0),'k.',fD,db(GD0),'k')
    ylabel('Amplitude (dB)')
    xlabel('f (Hz)')
    title(strcat('M=',num2str(MAll(k))))
    axis([0 60 -60 1])
    DG_SetFontSize(12)
    DG_SetTraceTint(30,2,FigNum,k)
    DG_SetTraceWidth(2,'*',FigNum,k)
end

% export the results
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation313.pdf', gcf); 
 
% print -dpdf FigExcitation313
 
 