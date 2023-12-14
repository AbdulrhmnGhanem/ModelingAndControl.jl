% Chapter 3 Exercise 38
% Revealing the nature of the leakage effects in FRF measurements
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all

% define the input variables
M=20;                       % number of blocks (sub-records)
NPeriodAll=[125 250 500]    % length of a block, program written for 3 lengths
N=max(NPeriodAll)*M;        % Number of data points in the measurement
fs=min(NPeriodAll);         % sampling frequency s.t. smallest period = 1 s
Ampl=1                      % RMS-value of the excitation
NTrans=1024;                % transient to eliminate initial effects at beginning simulation
t=[0:N-1]/fs;               % time axis

[b,a]=cheby1(2,10,0.2);     % test system

% generate excitation signal   
u=randn(N+NTrans,1)*Ampl;   % white noise excitation
    
% system response    
y=filter(b,a,u);                 % response of the system
y(1:NTrans)=[];u(1:NTrans)=[];   % eliminate transients of the simulation

% processing for different block lenghts
for r=length(NPeriodAll):-1:1
NPeriod=NPeriodAll(r);      % select the actual block length

Lines=[1:NPeriod/2+1];      % select the lines to be plotted
f=(Lines-1)/NPeriod*fs;     % frequency vector
M=round(N/NPeriod);         % number of blocks
MAll(r)=M;
          
    u=reshape(u,NPeriod,M); % put each block in a collumn
    U=fft(u);U=U(Lines,:);  % FFT-analysis
    
    G=freqz(b,a,2*pi*f/fs); % FRF
    
    YMod=zeros(NPeriod,M);  % initialization
    for s=1:M        % calculate response of block via FRF
        YMod(Lines,s)=G(:).*U(:,s);         % modelled output, assuming periodic repition
        yMod(:,s)=2*real(ifft(YMod(:,s)));
    end
    
    yModAll(:,r)=yMod(:);     % concatenate to a single vector
    eModAll(:,r)=yMod(:)-y(:);% error of the periodically reconstructed output
    clear yMod
    
end


% plot the results
t=[0:N-1];                     % time vector 

FigNum=10
figure(FigNum)
clf
n=3;                            % number of block lengths that are analysed

subplot(2,3,2)                  % plot the time signal
plot(t,y,'k')
ylabel('Output')
DG_SetFontSize(12)
for k=1:n                       % plot the transient errors
    subplot(2,3,k+3)
    tAll=reshape(t,NPeriodAll(k),MAll(k));
    e=eModAll(:,k);eAll=e(tAll+1);
    plot(tAll,eAll,'k')
    axis([0 4000 -1 1 ])
    if k==1; ylabel('Error');end
    DG_SetFontSize(12)
    for j=2:2:MAll(k)
    DG_SetTraceTint(30,j,FigNum,k+1)
    end
    DG_SetTraceWidth(2,'*',FigNum,k+1)
end
subplot(2,3,5)
xlabel('Sample number')
     
% Export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitation310.pdf', gcf);
 