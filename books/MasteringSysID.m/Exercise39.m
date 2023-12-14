% Chapter 3 Exercise 39
% FRF-measurement using a noise excitation and a Hanning window
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
% Compare rectangular and Hanning window for number of averages
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
y=reshape(y,NPeriod,M);          % put each period in a collumn
u=reshape(u,NPeriod,M);
    
% hanning window
w=hanning(NPeriod,'periodic');    % the Hanning window
wAll=kron(w(:),ones(1,M));        % repeated window to process all sub-records at once

Y=fft(y)/sqrt(NPeriod/2);Y=Y(Lines,:);                     % spectral analysis, block per block
U=fft(u)/sqrt(NPeriod/2);U=U(Lines,:);                     %        using rectangular window

YHann=fft(y.*wAll)/sqrt(NPeriod/2);YHann=YHann(Lines,:);   % spectral analysis, block per block
UHann=fft(u.*wAll)/sqrt(NPeriod/2);UHann=UHann(Lines,:);   % using Hanning window
    
for k=length(MAll):-1:1           % loop over number of blocks
    UU(:,k)=mean(abs(U(:,1:MAll(k))).^2,2);                 % estimate SUU rectangular
    YU(:,k)=mean(Y(:,1:MAll(k)).*conj(U(:,1:MAll(k))),2);   % estimate SYU
    G(:,k)=YU(:,k)./UU(:,k);                                % FRF rectangular

    UUHann(:,k)=mean(abs(UHann(:,1:MAll(k))).^2,2);         % estimate SUU Hanning
    YUHann(:,k)=mean(YHann(:,1:MAll(k)).*conj(UHann(:,1:MAll(k))),2);
                                                            % estimate SYU
    GHann(:,k)=YUHann(:,k)./UUHann(:,k);                    % FRF Hanning
end
    
% plot the results: FRF for increasing number of averaged blocks
FigNum=1 
figure(FigNum)   % using the Hanning window

G0=freqz(b,a,f*2*pi/fs);
clf
subplot(2,3,1)
plot(f,db(G0),'k',f,db(G(:,1)-G0),'.k',f,db(GHann(:,1)-G0),'.k'),title('M=1')
ylabel('Amplitude (dB)')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,3,2)
plot(f,db(G0),'k',f,db(G(:,2)-G0),'.k',f,db(GHann(:,2)-G0),'.k'),title('M=4')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,3,3)
plot(f,db(G0),'k',f,db(G(:,3)-G0),'.k',f,db(GHann(:,3)-G0),'.k'),title('M=16')

axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,3,4)

plot(f,db(G0),'k',f,db(G(:,4)-G0),'.k',f,db(GHann(:,4)-G0),'.k'),title('M=256')
axis([0 60 -80 0])
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
DG_SetTraceTint(30,2,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

subplot(2,3,5)

plot(f,db(G0),'k',f,db(G(:,5)-G0),'.k',f,db(GHann(:,5)-G0),'.k'),title('M=1024')
xlabel('f (Hz)')

axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,5)
DG_SetTraceWidth(2,'*',FigNum,5)
DG_SetFontSize(12)

subplot(2,3,6)

plot(f,db(G0),'k',f,db(G(:,6)-G0),'.k',f,db(GHann(:,6)-G0),'.k'),title('M=4096')
xlabel('f (Hz)')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)
DG_SetFontSize(12)

%DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
%DG_MakePDF('FigExcitation307Hann1.pdf', gcf); 
    
    
%
%%
%%%
% Part 2
%
% Compare rectangular and Hanning window for different window lengths, and fixed
% number of averages
%%%
%%
%


% define the input variables
M=128;  % number of blocks to be processed
NPeriodAll=[128 256 512 1024 2048 4096];       % length of a block
clear f UU YU G0 G UUHann YUHann GHann


for k=length(NPeriodAll):-1:1
NPeriod=NPeriodAll(k);  % select the actual block length    
N=NPeriod*M;       % Number of data points in the measurement
fs=128         % sampling frequency is 128 Hz
Ampl=1             % RMS-value of the multi sine
fGen=0.3*fs        % bandwidth excitation
NTrans=1024    % to eliminate transient effects

t=[0:NPeriod-1]/fs;   % time vector 1 period
tAll=[0:N-1]/fs;      % time vector full sequence
Lines=[1:NPeriod/2];% select the lines to be plotted
f{k}=[Lines-1]'/NPeriod*fs; % frequency vector



% generate excitation signal
    
u=randn(N+NTrans,1)*Ampl;
u=filter(bGen,aGen,u);

% system response    

y=filter(b,a,u);                 % response of the system
y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients

% hanning window
w=hanning(NPeriod,'periodic');   % the Hanning window
wAll=kron(w(:),ones(1,M));   % repeated window to process all sub-records at once

% processing: steady state and transient    
y=reshape(y,NPeriod,M);             % put each period in a collumn
u=reshape(u,NPeriod,M);

Y=fft(y)/sqrt(NPeriod/2);Y=Y(Lines,:);            % spectral analysis, block per block
U=fft(u)/sqrt(NPeriod/2);U=U(Lines,:);            %        using rectangular window

YHann=fft(y.*wAll)/sqrt(NPeriod/2);YHann=YHann(Lines,:);   % spectral analysis, block per block
UHann=fft(u.*wAll)/sqrt(NPeriod/2);UHann=UHann(Lines,:);   % using Hanning window

% process the data
UU{k}=mean(abs(U(:,1:M)).^2,2);
YU{k}=mean(Y(:,1:M).*conj(U(:,1:M)),2);
G{k}=YU{k}./UU{k};

UUHann{k}=mean(abs(UHann(:,1:M)).^2,2);
YUHann{k}=mean(YHann(:,1:M).*conj(UHann(:,1:M)),2);
GHann{k}=YUHann{k}./UUHann{k};


% true system G0

G0{k}=freqz(b,a,f{k}*2*pi/fs);
end



% plot the results: FRF for changing block length, fixed number of averages
FigNum=2 
figure(FigNum)   % using the Hanning window


clf
subplot(2,3,1)
plot(f{1},db(G0{1}),'k',f{1},db(G{1}-G0{1}),'.k',f{1},db(GHann{1}-G0{1}),'.k'),title('N=128')
ylabel('Amplitude (dB)')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,3,2)
plot(f{2},db(G0{2}),'k',f{2},db(G{2}-G0{2}),'.k',f{2},db(GHann{2}-G0{2}),'.k'),title('N=256')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,2)
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,3,3)
plot(f{3},db(G0{3}),'k',f{3},db(G{3}-G0{3}),'.k',f{3},db(GHann{3}-G0{3}),'.k'),title('N=512')

axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,3,4)

plot(f{4},db(G0{4}),'k',f{4},db(G{4}-G0{4}),'.k',f{4},db(GHann{4}-G0{4}),'.k'),title('N=1024')
axis([0 60 -80 0])
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
DG_SetTraceTint(30,2,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

subplot(2,3,5)

plot(f{5},db(G0{5}),'k',f{5},db(G{5}-G0{5}),'.k',f{5},db(GHann{5}-G0{5}),'.k'),title('N=2048')
xlabel('f (Hz)')
ylabel('Amplitude (dB)')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,5)
DG_SetTraceWidth(2,'*',FigNum,5)
DG_SetFontSize(12)

subplot(2,3,6)

plot(f{6},db(G0{6}),'k',f{6},db(G{6}-G0{6}),'.k',f{6},db(GHann{6}-G0{6}),'.k'),title('N=4096')
xlabel('f (Hz)')
axis([0 60 -80 0])
DG_SetTraceTint(30,2,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)
DG_SetFontSize(12)

DG_SetTraceTint(30,2,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)
DG_SetFontSize(12)
 
% Export the results
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation307Hann2.pdf', gcf); 