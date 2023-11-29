% Chapter 3 Exercise 36
% Study of a multisine response of a linear system: transients and steady
% state
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables

NPeriod=128             % period length multisine
M=3                     % generate 3 periods
N=NPeriod*M;            % Number of data points in the measurement
fs=NPeriod*2            % sampling frequency s.t. 1 period = 1 s
Ampl=1                  % RMS-value of the multi sine


t=[0:NPeriod-1]/fs;     % time vector 1 period
tAll=[0:N-1]/fs;        % time vector full sequence
Lines=1+[1: floor(0.8*NPeriod/2)]';% select the lines to be plotted
f=(Lines-1)/NPeriod*fs; % frequency vector

[b,a]=cheby1(2,10,0.2); % test system

% generate full band multisine
ExcitationLines=1+ [1: floor(0.8*NPeriod/2)];    % excite till 80% full band
U=zeros(NPeriod,1);     % prepare calculation multisine via IFFT
U(ExcitationLines)=exp(j*2*pi*rand(length(ExcitationLines),1));  % random phase multisine Fourier coeff.
u=real(ifft(U));u=u/max(abs(u))*Ampl;                             
                        % generate 1 period with normalized peak value
u=kron(ones(M,1),u);    % repeat M times 
    
% system response    
yM=filter(b,a,u);       % response of the system
    
% processing: extract steady state and transient    
y=reshape(yM,NPeriod,M);                % put each period in a collumn
u=reshape(u,NPeriod,M);
ySteadyState=kron(ones(1,M),y(:,3));    % repeat the last period fSine times
yTransient=y-ySteadyState;              % calculate the transient

% Plot the results
FigNum=3;
figure(FigNum)
clf
subplot(1,3,1)          % plot time signals
i1=[1:NPeriod+1];
i2=[NPeriod+1:2*NPeriod+1];
i3=[2*NPeriod+1:3*NPeriod];
plot(tAll(i1),y(i1),'k',tAll(i2),y(i2),'k',tAll(i3),y(i3),'k'),ylabel('Output'),xlabel('Time (s)')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(2,'*',FigNum,1)
DG_SetFontSize(12)

subplot(1,3,2)          % plot the transient
plot(tAll,yTransient(:),'k'),ylabel('Transient'),xlabel('Time (s)')
DG_SetTraceWidth(2,'*',FigNum,2)
DG_SetFontSize(12)

subplot(1,3,3)          % plot the transient on a logarithmic scale
semilogy(tAll,abs(yTransient(:)),'k'),ylabel('abs(Transient)'),xlabel('Time (s)')
DG_SetTraceWidth(2,'*',FigNum,3)
DG_SetFontSize(12)
axis([0 1.5 1e-10 1])
    
% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation303.pdf', gcf);
  
    
% extracting the FRF
U=fft(u)/sqrt(NPeriod);     % FFT input
Y=fft(y)/sqrt(NPeriod);     % FFT output
G=Y(Lines,:)./U(Lines,:);   % calculate FRF

% Plot the results
FigNum=4;	
figure(FigNum)

subplot(2,2,1)                  % input spectrum
plot(f,abs(U(Lines,end)),'.k'),title('Input (a)')
ylabel('Amplitude')
axis([0 NPeriod 0 1])
DG_SetFontSize(12)

subplot(2,2,2)
G0=freqz(b,a,2*pi*f/fs);         % FRF, G0, error  1st period
plot(f,db(G0),'k',f,db(G(:,1)),'.k',f,db(G(:,1)-G0),'k'),title('FRF (b)')
DG_SetTraceTint(30,3,FigNum,2)
DG_SetTraceWidth(2,3,FigNum,2) 
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
axis([0 NPeriod -60 0])
DG_SetFontSize(12)

subplot(2,2,3)
G0=freqz(b,a,2*pi*f/fs);         % FRF, G0, error  2nd period
plot(f,db(G(:,2)),'.k',f,db(G(:,2)-G0),'k'),title('FRF (c)')
DG_SetTraceTint(30,2,FigNum,3)
DG_SetTraceWidth(2,2,FigNum,3)
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
axis([0 NPeriod -200 0])
DG_SetFontSize(12)

 subplot(2,2,4)
G0=freqz(b,a,2*pi*f/fs);         % FRF, G0, error  3th period
plot(f,db(G(:,3)),'.k',f,db(G(:,3)-G0),'k'),title('FRF (d)')
DG_SetTraceTint(30,2,FigNum,4)
DG_SetTraceWidth(2,2,FigNum,4)
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
axis([0 NPeriod -200 0])
    
% Export the plot   
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation304.pdf', gcf);

