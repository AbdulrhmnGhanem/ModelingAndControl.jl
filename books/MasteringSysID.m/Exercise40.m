    
% Chapter 3 Exercise 40
% FRF-measurement using a noise excitation and a diff window
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

%
%%
%%%
% Part 1
%
% Compare rectangular and diff window for number of averages
% and a fixed block length
%
%%%
%%
%

clear all

% define the input variables
MAll=[1 4 16 256 1024 4096];    % number of blocks to be processed
M=max(MAll);
NPeriod=1024                    % length of a block
N=NPeriod*M;                    % Number of data points in the measurement
fs=128                          % sampling frequency
Ampl=1                          % RMS-value of the multi sine
fGen=0.3*fs                     % bandwidth excitation
NTrans=1024                     % to eliminate transient effects



t=[0:NPeriod-1]/fs;             % time vector of 1 period
tAll=[0:N-1]/fs;                % time vector of full sequence
Lines=[1:NPeriod/2];            % select the lines to be plotted
f=[Lines-1]'/NPeriod*fs;        % frequency vector

[b,a]=cheby1(2,10,0.2);         % test system
[bGen,aGen]=butter(2,fGen/fs*2);% generator filter to generate colored random excitation

% generate excitation signal   
u=randn(N+NTrans,1)*Ampl;       % generate white noise sequence
u=filter(bGen,aGen,u);          % colored input noise
    
% system response    
y=filter(b,a,u);                 % response of the system
y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients

% processing: steady state and transient    
y=reshape(y,NPeriod,M);                         % put each period in a collumn
u=reshape(u,NPeriod,M);

Lines1=[1:NPeriod/2+1];                        % 1 additional line selected to apply diff window
Y=fft(y)/sqrt(NPeriod/2);Y=Y(Lines1,:);         % FFT
Y=diff(Y);                                      % apply diff window
U=fft(u)/sqrt(NPeriod/2);U=U(Lines1,:);         % FFT
U=diff(U);                                  % apply diff window
f=f+0.5/NPeriod*fs;                             % frequency compensation diff-window


    


for k=length(MAll):-1:1           % loop over number of blocks
    UU(:,k)=mean(abs(U(:,1:MAll(k))).^2,2);         % estimate SUU
    YU(:,k)=mean(Y(:,1:MAll(k)).*conj(U(:,1:MAll(k))),2);   % esimtat SYU
    G(:,k)=YU(:,k)./UU(:,k);
end
    
% Plot the results: FRF for increasing number of averaged blocks
FigNum=7 
figure(FigNum)   % using the diff window
clf

G0=freqz(b,a,f*2*pi/fs);
clf
subplot(2,3,1)
plot(f,db(G(:,1)),'k',f,db(G0),'k',f,db(G(:,1)-G0),'.k'),title('M=1')
ylabel('Amplitude (dB)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,3,2)
plot(f,db(G(:,2)),'k',f,db(G0),'k',f,db(G(:,2)-G0),'.k'),title('M=4')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,3,3)
plot(f,db(G(:,3)),'k',f,db(G0),'k',f,db(G(:,3)-G0),'.k'),title('M=16')
ylabel('Amplitude (dB)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,3,4)
xlabel('f (Hz)')
plot(f,db(G(:,4)),'k',f,db(G0),'k',f,db(G(:,4)-G0),'.k'),title('M=256')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

subplot(2,3,5)
xlabel('f (Hz)')
plot(f,db(G(:,5)),'k',f,db(G0),'k',f,db(G(:,5)-G0),'.k'),title('M=1024')
xlabel('f (Hz)')
ylabel('Amplitude (dB)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,5)
DG_SetTraceWidth(2,'*',FigNum,5)
DG_SetFontSize(12)

subplot(2,3,6)
xlabel('f (Hz)')
plot(f,db(G(:,6)),'k',f,db(G0),'k',f,db(G(:,6)-G0),'.k'),title('M=4096')
xlabel('f (Hz)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)
DG_SetFontSize(12)
 
% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation307.pdf', gcf);   





%
%%
%%%
% Part 2
%
% Compare rectangular and diff window for different window lengths, and fixed
% number of averages
%%%
%%
%



% define the input variables
NPeriodAll=[128 128*4 128*16 128*64];  % number of blocks to be processed
M=128;
N=NPeriodAll(end)*M;            % Number of data points in the measurement


% generate excitation signal   
u=randn(N+NTrans,1)*Ampl;       % generate white noise sequence
u=filter(bGen,aGen,u);          % colored input noise
    
% system response    

    y=filter(b,a,u);                 % response of the system
    y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients
    
    
 % process the data
 
 for r=length(NPeriodAll):-1:1
    NPeriod=NPeriodAll(r)
    M1=N/NPeriod;
    u=reshape(u,NPeriod,M1);            % 1 block per collumn
    y=reshape(y,NPeriod,M1);
    
    Lines=[1:NPeriod/2];                % 1 additional for the diff window
    Y=fft(y)/(NPeriod/2);YR=Y(Lines,:); % spectral analysis
    U=fft(u)/(NPeriod/2);UR=U(Lines,:);
    Y=diff(Y);U=diff(U);
    YD=Y(Lines,:);UD=U(Lines,:);
    
    fR{r}=(Lines-1)'/NPeriod;           % frequency vector rectangular window
    fD{r}=fR{r}+0.5/NPeriod;            % frequency vector diff window
    
    % process with rectangular window
    UU=mean(abs(UR(:,1:M)).^2,2);       % SUU  
    YU=mean(YR(:,1:M).*conj(UR(:,1:M)),2); %SYU
    GR{r}=YU./UU;                       % FRF

    % process with Diff window
    UU=mean(abs(UD(:,1:M)).^2,2);       % SUU
    YU=mean(YD(:,1:M).*conj(UD(:,1:M)),2); % SYU
    GD{r}=YU./UU;                       % SYU
        

    G0R{r}=freqz(b,a,fR{r}*2*pi);       % Store the exact G0
    G0D{r}=freqz(b,a,fD{r}*2*pi);
      
    
end   % loop over block lengths


% Plot the results: FRF for changing block length, fixed number of averages
FigNum=8;
figure(FigNum)
clf
subplot(2,2,1)
plot(fR{1},db(G0R{1}),'k',fR{1},...
          db(GR{1}-G0R{1}),'.k',fD{1},db(GD{1}-G0D{1}),'.k'),title('N=128')
ylabel('Amplitude (dB)')
axis([0 0.5 -100 0])
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,2,2)
plot(fR{2},db(G0R{2}),'k',...
         fR{2},db(GR{2}-G0R{2}),'.k',fD{2},db(GD{2}-G0D{2}),'.k'),title('N=512')
axis([0 0.5 -100 0])
DG_SetTraceTint(30,2,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,2,3)
plot(fR{3},db(G0R{3}),'k',...
          fR{3},db(GR{3}-G0R{3}),'.k',fD{3},db(GD{3}-G0D{3}),'.k'),title('N=2048')
ylabel('Amplitude (dB)')
xlabel('f/fs')
axis([0 0.5 -100 0])
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,2,4)
plot(fR{4},db(G0R{4}),'k',...
           fR{4},db(GR{4}-G0R{4}),'.k',fD{4},db(GD{4}-G0D{4}),'.k'),title('N=8192')
xlabel('f/fs')
axis([0 0.5 -100 0])
DG_SetTraceTint(30,2,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation308.pdf', gcf); 
