% Chapter 3 Exercise 37
% FRF-measurement using a noise excitation and a rectangular window
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables
MAll=[1 4 16 256 1024 4096];    % number of blocks to be processed
M=max(MAll);
NPeriod=128                     % length of a block
N=NPeriod*M;                    % Number of data points in the measurement
fs=NPeriod                      % sampling frequency s.t. 1 period = 1 s
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
Y=fft(y)/sqrt(NPeriod/2);Y=Y(Lines,:);          % spectral analysis, block per block
U=fft(u)/sqrt(NPeriod/2);U=U(Lines,:);
    
for k=length(MAll):-1:1                      % loop over period numbers
    UU(:,k)=mean(abs(U(:,1:MAll(k))).^2,2);                % estimate SUU
    YU(:,k)=mean(Y(:,1:MAll(k)).*conj(U(:,1:MAll(k))),2);  % estimate SYU
    G(:,k)=YU(:,k)./UU(:,k);                               % estimate FRF
end
    
% Plot the results: the FRF for increasing number of averaged blocks   
FigNum=5;
figure(FigNum)
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
ylabel('Amplitude (dB)')
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
plot(f,db(G(:,4)),'k',f,db(G0),'k',f,db(G(:,4)-G0),'.k'),title('M=256')
ylabel('Amplitude (dB)'),xlabel('f (Hz)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,4)
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

subplot(2,3,5)
plot(f,db(G(:,5)),'k',f,db(G0),'k',f,db(G(:,5)-G0),'.k'),title('M=1024')
ylabel('Amplitude (dB)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,5)
DG_SetTraceWidth(2,'*',FigNum,5)
DG_SetFontSize(12)
xlabel('f (Hz)')

subplot(2,3,6)
plot(f,db(G(:,6)),'k',f,db(G0),'k',f,db(G(:,6)-G0),'.k'),title('M=4096')
ylabel('Amplitude (dB)')
axis([0 60 -60 0])
DG_SetTraceTint(30,1,FigNum,6)
DG_SetTraceWidth(2,'*',FigNum,6)
DG_SetFontSize(12)
xlabel('f (Hz)')

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation305.pdf', gcf);
     
 %   print -dpdf FigExcitation305
    
FigNum=6;
figure(FigNum) % Evolution sample power spectrum, impact number of averages

clf
subplot(2,2,1)
plot(f,db(UU(:,1)/2),'k'),title('M=1')
ylabel('Amplitude (dB)')
axis([0 60 -40 20])
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(2,2,2)
plot(f,db(UU(:,2)/2),'k'),title('M=4')
axis([0 60 -40 20])
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(2,2,3)
plot(f,db(UU(:,3)/2),'k'),title('M=16')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')
axis([0 60 -40 20])
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)

subplot(2,2,4)
xlabel('Frequency (Hz)')
plot(f,db(UU(:,4)/2),'k'),title('M=256')
xlabel('f (Hz)')
axis([0 60 -40 20])
DG_SetTraceWidth(2,'*',FigNum,4)
DG_SetFontSize(12)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation306.pdf', gcf);

 
 