% Chapter 3 Exercise 53
% Measuring the FRM using multisine excitations
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 29 November 2010

clear all

% define the input variables
NPeriod=64          % period length multisine
M=2                 % generate 2 periods
N=NPeriod*M;        % Number of data points in the measurement
fs=NPeriod*2        % sampling frequency
Ampl=1              % RMS-value of the multi sine
fMax=0.4            % max freq of the excitation

t=[0:NPeriod-1]/fs;                 % time vector 1 period
Lines=1+[1: floor(0.8*NPeriod/2)]'; % select the lines to be plotted
f=(Lines-1)/NPeriod*fs;             % frequency vector

% define the MIMO system
[b{1,1} a{1,1}]=cheby1(2,10,0.2);   % G11
[b{1,2} a{1,2}]=cheby1(2,20,0.5);   % G12
[b{2,1} a{2,1}]=butter(4,0.6);      % G21
[b{2,2} a{2,2}]=cheby1(2,15,0.8);   % G22

% calculate FRM G0  
for r=1:2   
    for s=1:2
        G0(:,r,s)=freqz(b{r,s},a{r,s},f*2*pi/fs);
    end
end

% generate full band multisine for input u      
ExcitationLines=1+ [1: floor(0.4*NPeriod)];                     % excite till 80% full band
U=zeros(NPeriod,1);                                             % initialize U
U(ExcitationLines)=exp(j*2*pi*rand(length(ExcitationLines),1)); % random phase multisine
u0=real(ifft(U));u0=u0/max(abs(u0))*Ampl;                       % generate 1 period with normalized peak value
u=kron(ones(M,1),u0);                                           % repeat M times 

% Gnerate two experiments to identify a 2x2 mimo system    
% Experiment 1   [u;u]
for r=1:2
   y{r,1}=filter(b{r,1},a{r,1},u)+filter(b{r,2},a{r,2},u);      % Gr1 u + Gr2 u
end

% Experiment 2    [u;-u]
for r=1:2
   y{r,2}=filter(b{r,1},a{r,1},u)+filter(b{r,2},a{r,2},-u);     % Gr1 u + Gr2 (-u)
end
       
% processing: take out the last period M as steady state solution  

U=fft(u0)/sqrt(NPeriod);    % process the input
U=U(Lines);                 % select the input lines

for r=1:2    % process the outputs
    for s=1:2
	y{r,s}=reshape(y{r,s},NPeriod,M);x=y{r,s};y{r,s}=x(:,M);    % take out the last period
    X=fft(y{r,s})/sqrt(NPeriod);                                % FFT analysis
    Y{r,s}=X(Lines);                                            % select the Lines that will be used
    end
end


for k=length(Lines):-1:1    % estimate the FRM at frequency k
  Yk=[Y{1,1}(k) Y{1,2}(k);Y{2,1}(k) Y{2,2}(k)];
  Uk=[U(k) U(k);U(k) -U(k)];
  GEst(k,:,:)=Yk/Uk;
end

% plot the results
FigNum=1;
figure(FigNum)
clf

counter=0;
for r=1:2
    for s=1:2
        counter=counter+1;
        subplot(2,2,counter)
        plot(f,db(G0(:,r,s)),'k',f,db(G0(:,r,s)-GEst(:,r,s)),'.k')
        if s==1;ylabel('Amplitude (dB)');end
        if r==2;xlabel('f (Hz)');end

    end
end

DG_SetTraceTint(30,2,FigNum,'*')
DG_SetTraceWidth(2,'*',FigNum,'*')
DG_SetFontSize(12) 

% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitationMIMO1.pdf', gcf);
