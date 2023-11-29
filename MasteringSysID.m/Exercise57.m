% Chapter 3 Exercise 57
% Measuring the FRM using noise excitations and a Hanning window
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 30 November 2010

clear all

% define the input variables
NRep=100            % number of repeated simulations
NBlock=512          % Block length of the sub-records in the processing
M=50                % number of sub-records
NTrans=1024         % eliminate transient effects simulation
fs=1                % sampling frequency = 1 Hz
Ampl=1              % RMS-value of the multi sine
fMax=0.4            % max freq of the excitation
nu=2                % number of inputs is froozen to two in this program

% define the MIMO system
[b{1,1} a{1,1}]=cheby1(2,10,0.2);       % G11
[b{1,2} a{1,2}]=cheby1(2,20,0.5);       % G12
[b{2,1} a{2,1}]=butter(4,0.6);          % G21
[b{2,2} a{2,2}]=cheby1(2,15,0.8);       % G22

% define the generation filter
[bGen aGen]=butter(4,0.8);             % generator filter

t=[0:NBlock-1]/fs;  % time vector 1 period
N=NBlock*M;         % Number of data points in the measurement
Lines=1+[1: floor(0.8*NBlock/2)]';      % select the lines FFT analysis
f=(Lines-1)/NBlock*fs;                  % frequency vector

% prepare the Hanning window
wHann=hann(NBlock,'periodic');     % Hanning window
wHann=kron(wHann,ones(1,M));       % matrix to be premultiplied with the 
                                   % reshaped signal vector for fast
                                   % calculations

% calculate FRM G0  
for r=1:2   
    for s=1:2
        G0(:,r,s)=freqz(b{r,s},a{r,s},f*2*pi/fs);
    end
end

for nRep=NRep:-1:1           % loop over the simulations
[nRep NRep]
% generate two filtered noise excitation signals u1 and u2       
u1=randn(N+NTrans,1);u1=filter(bGen,aGen,u1);
u2=randn(N+NTrans,1);u2=filter(bGen,aGen,u2);    
% system response   
y1=filter(b{1,1},a{1,1},u1)+filter(b{1,2},a{1,2},u2);       % G11 u1 + G12 u2
y2=filter(b{2,1},a{2,1},u1)+filter(b{2,2},a{2,2},u2);       % G21 u1 + G22 u2
    
% processing of the data 
u1(1:NTrans)=[]; u2(1:NTrans)=[]; %take out the steady state solution
y1(1:NTrans)=[]; y2(1:NTrans)=[]; 

u1=reshape(u1,NBlock,M);u2=reshape(u2,NBlock,M);    % break data in subrecords
y1=reshape(y1,NBlock,M);y2=reshape(y2,NBlock,M);

U1=fft(u1.*wHann)/sqrt(NBlock);U2=fft(u2.*wHann)/sqrt(NBlock);    % calculate windowed FFT
Y1=fft(y1.*wHann)/sqrt(NBlock);Y2=fft(y2.*wHann)/sqrt(NBlock);


U1=U1(Lines,:);U2=U2(Lines,:);            % select the lines of interest
Y1=Y1(Lines,:);Y2=Y2(Lines,:); 

for k=length(Lines):-1:1
    SUU(k,:,:)=[U1(k,:);U2(k,:)]*[U1(k,:);U2(k,:)]';       % calculate SUU
    SYU(k,:,:)=[Y1(k,:);Y2(k,:)]*[U1(k,:);U2(k,:)]';       % calculate SYU
    SYY(k,:,:)=[Y1(k,:);Y2(k,:)]*[Y1(k,:);Y2(k,:)]';       % calculate SYY
end


for k=length(Lines):-1:1    % estimate the FRM at frequency k
  GEst(k,nRep,:,:)=squeeze(SYU(k,:,:))/squeeze(SUU(k,:,:));
end


for k=length(Lines):-1:1    % estimate the cov(vec(Gest)) at frequency k
  SUUinv=inv(squeeze(SUU(k,:,:)));
  C=squeeze(SYY(k,:,:))-squeeze(SYU(k,:,:))*SUUinv*squeeze(SYU(k,:,:))';
  C=1/(M-nu)*kron(SUUinv,C);
  
  VarG(k,nRep,:,:)=[C(1,1) C(3,3); C(2,2) C(4,4)];   
                            % variance [G11 G12; G21 G22]
end  % loop over the frequencies

end   % loop over nRep

% process the data of the repeated simulations
for k=length(Lines):-1:1  % loop over the frequencies
  GEstm(k,:,:)=mean(GEst(k,:,:,:));     % mean value of G over the realizations
  VarGEst(k,:,:)=var(GEst(k,:,:,:));    % estimated variance G over the realizations
  VarGm(k,:,:)=mean(VarG(k,:,:,:));     % mean theoretic variance over the realizations
end % loop over the frequencies
  

% plot the figures  
    FigNum=1;
    figure(FigNum)
    clf
    
    counter=0;
    for r=1:2
        for s=1:2
            counter=counter+1;
            subplot(2,2,counter)
            plot(f,db(G0(:,r,s)),'k',...
                 f,db(G0(:,r,s)-GEstm(:,r,s)),'.k',f,db(VarGEst(:,r,s))/2,'k',...
                 f,db(VarGm(:,r,s))/2,'k')
             axis([0 0.4 -80 5])
            if s==1;ylabel('Amplitude (dB)');end
            if r==2;xlabel('f (Hz)');end
        end
    end
    
DG_SetTraceTint(30,2,FigNum,'*')
DG_SetTraceTint(60,3,FigNum,'*')
DG_SetTraceTint(1,4,FigNum,'*')
DG_SetTraceWidth(2,'*',FigNum,'*')
DG_SetTraceWidth(4,3,FigNum,'*')
DG_SetTraceWidth(1,4,FigNum,'*')
DG_SetFontSize(12) 
 
% export the plots
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigExcitationMIMO3.pdf', gcf);

