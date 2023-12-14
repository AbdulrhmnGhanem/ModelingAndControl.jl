% Chapter 1 Exercise 12
% Noise on input and output: the instrumental variables method applied on the resistor estimate 
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

clear all
N=5000           % number of datapoints
NRep=10000       % number of repeated experiments
R0=1000          % exact value resistor
iMax=0.01        % maximum current
fGen=0.1         % cutt-off freq. filter for current generator (1=fs/2)
CuttOffAll=[0.999 0.95 0.6];    % cutt off freq. noise filter  1=fs/2
NTrans=1000;     % to eliminate transient effects in the filter operations

% Generate the exact current
[bGen,aGen]=butter(1,fGen);                           % current generating filter
i0=randn(N+NTrans,1);i0=filter(bGen,aGen,i0);         % generate current
i0(1:NTrans)=[];i0=i0*iMax/std(i0);                   % scale the filtered signal

GenCorr=xcorr(i0,11,'unbiased');    % autocorrelation current --> used for plotting only!
GenCorr(1:11)=[];                   % only the positive lags are used


u0=R0*i0;   % exact voltage

% Part1:  changing the bandwidth of the noise on the current, 
%         fixed shift of 1 sample


nShift=1;   % shift IV over nShift points
clf
for r=length(CuttOffAll):-1:1                           % loop over the different noise filters
    
    [b,a]=butter(2,CuttOffAll(r));                      % create the disturbing noise filter
    G(:,r)=abs(freqz(b,a,N,'whole'));                   % used to plot the figures
    
    for s=1:NRep                                        % loop over the different realizations
        iNoise=randn(N+NTrans,1);
        iNoise=filter(b,a,iNoise);iNoise(1:NTrans)=[];  % filtered disturbing noise on the current
        iNoise=iNoise/std(iNoise)*iMax;                 % scaling
        
        uNoise=randn(N,1);                             % white disturbing noise on the voltage
        
        i=i0+iNoise;                                   % add noise to the current measurements
        u=u0+uNoise;                                   % add noise to the voltage measurements
        
        RLS(s,r)=(i'*u)/(i'*i);                        % Least squares estimate
        
        iShift=i;iShift(1:nShift)=[];                  % generate shifted vectors
        u(end+1-nShift:end)=[];i(end+1-nShift:end)=[];
        
        RIV(s,r)=(iShift'*u)/(iShift'*i);              %  calculated IV
        
    end
    
    z=xcorr(iNoise,11,'unbiased');nCorr(:,r)=z(12:end); % estimate autocorrelation noise
end


% Plot the first series of results
n=length(CuttOffAll);
FigNum=1;
figure(FigNum)

f=[0:(N-1)/2]'/N;
BinWidth=1;X=[0:BinWidth:1500];   % set the parameters of the HIST call to estimate the pdf

for r=1:n    % loop over the different noise filters
    pdfRLS(:,r)=hist(RLS(:,r),X)/BinWidth/NRep;     % estimated pdf
    pdfIV(:,r)=hist(RIV(:,r),X)/BinWidth/NRep;      % estimated pdf
end

figure(FigNum),clf
    plot(X,pdfRLS,'k',X,pdfIV,'k')
    axis([0 1500 0 0.1])
    
    DG_SetFontSize(8)
    ylabel('pdf(R)')
    DG_SetTraceTint(30,5,FigNum,1)
    DG_SetTraceWidth(1,'*',FigNum,1)

% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIVnr1a.pdf', gcf);

% plot the correlation functions of i0 and noise on i0
FigNum=2
figure(FigNum),clf
GGen=abs(freqz(bGen,aGen,N,'whole'));                 % generator filter characteristic

for r=1:n    
    subplot(2,n,r)                                % plot the correlation functions 
    plot([0:10],real(GenCorr(1:11)),'+-k',[0:10],real(nCorr(1:11,r)),'+-k'),hold on
    plot(nShift,real(GenCorr(nShift+1)),'ok',nShift,real(nCorr(nShift+1,r)),'ok'),hold off
    DG_SetFontSize(12)
    if r==1; ylabel('Auto-correlation'); end
    xlabel('Lag number')
    axis([0 10 -0.2e-4 1e-4])
    
    subplot(2,n,r+n)                               % plot the filter characteristics
    plot(f,db(GGen(1:N/2)),'-k',f,db(G(1:N/2,r)),'-k')
    DG_SetFontSize(12)
    axis([0 0.5 -30 10])
    if r==1; ylabel('Filter (dB)'); end
    xlabel('f/fs')
end

for r=1:n   % set gray scales of last row of plots
    DG_SetTraceTint(30,2,FigNum,r)
    DG_SetTraceWidth(1,'*',FigNum,r)
end

for r=1:n   % set gray scales of last row of plots
    DG_SetTraceTint(30,2,FigNum,r+3)
    DG_SetTraceWidth(1,'*',FigNum,r+3)
end

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIVnr1b.pdf', gcf);
    
%print -dpdf FigIVnr1



% Part2:  constant bandwidth of the noise on the current, 
%         changing the shift
clear all
N=5000           % number of datapoints
fGen=0.1         % cutt-off freq. current generator (2=fs)
NRep=10000       % number of repeated experiments
R0=1000          % exact value resistance
iMax=0.01        % maximum current
NTrans=1000;     % to eliminate transient effects in the filter operations
CuttOffAll=0.6   % cut off frequency for the noise filter (2=fs)


nShiftAll=[1 2 5];   % shift IV over nShift points

[b,a]=butter(2,CuttOffAll(end)); % disturbing current noise filter: lowest cut-off freq.
G=abs(freqz(b,a,N,'whole'));     % noise filter FRF for plot

[bGen,aGen]=butter(1,fGen);                           % current generating filter
GGen=abs(freqz(bGen,aGen,N,'whole'));                 % filter characteristic for plot
i0=randn(N+NTrans,1);i0=filter(bGen,aGen,i0);         % generate undisturbed current i0
i0(1:NTrans)=[];i0=i0*iMax/std(i0);                   % scale the filtered signal
u0=R0*i0;                                             % exact voltage
GenCorr=xcorr(i0,11,'unbiased'); GenCorr(1:11)=[];    % autocorrelation inputv --> used for plotting only!

for r=length(nShiftAll):-1:1          % loop over the different shift values
    nShift=nShiftAll(r)
 
    for s=1:NRep                      % loop over the different realization
        iNoise=randn(N+NTrans,1);
        iNoise=filter(b,a,iNoise);iNoise(1:NTrans)=[]; % filtered disturbing current noise
        iNoise=iNoise/std(iNoise)*iMax;                % normalize RMS-value
        
        uNoise=randn(N,1);            % white disturbing noise on the voltage
        
        i=i0+iNoise;                  % add noise
        u=u0+uNoise;
        
        RLS(s,r)=(i'*u)/(i'*i);       % Least squares estimate
        
        iShift=i;iShift(1:nShift)=[];
        u(end+1-nShift:end)=[];i(end+1-nShift:end)=[]; 
        RIV(s,r)=(iShift'*u)/(iShift'*i);  % calculate IV
        
    end
     z=xcorr(iNoise,11,'unbiased');nCorr(:,r)=z(12:end);   % estimate autocorrelation noise for plot only
    
end
        
n=length(nShiftAll);





f=[0:(N-1)/2]'/N;
BinWidth=1;X=[0:BinWidth:1500];   % set HIST parameters to estimate pdf results

clear pdfRLS pdfIV
for r=1:n    % loop over the different noise filters
    pdfRLS(:,r)=hist(RLS(:,r),X)/BinWidth/NRep;   % estimated pdf
    pdfIV(:,r)=hist(RIV(:,r),X)/BinWidth/NRep;   % estimated pdf
end

% plot the results
FigNum=3;
figure(FigNum),clf     

figure(FigNum),clf
plot(X,pdfRLS,'k',X,pdfIV,'k')
axis([0 1500 0 0.05])
    
DG_SetFontSize(8)
ylabel('pdf(R)')
DG_SetTraceTint(30,5,FigNum,1)
DG_SetTraceWidth(1,'*',FigNum,1)
    
% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIVnr2a.pdf', gcf);

FigNum=FigNum+1;
figure(FigNum),clf
for r=1:n    
    subplot(3,n,r)                                % plot the correlation functions 
    nShift=nShiftAll(r);
    plot([0:10],real(GenCorr(1:11)),'+-k',[0:10],real(nCorr(1:11,r)),'+-k'),hold on
    plot(nShift,real(GenCorr(nShift+1)),'ok',nShift,real(nCorr(nShift+1,r)),'ok'),hold off
    DG_SetFontSize(12)
    if r==1; ylabel('Auto-correlation'); end
    xlabel('Lag number')
    axis([0 10 -0.2e-4 1e-4])
    
    subplot(3,n,r+n)                               % plot the filter characteristics
    plot(f,db(GGen(1:N/2)),'-k',f,db(G(1:N/2)),'-k')
    DG_SetFontSize(12)
    axis([0 0.5 -30 10])
    if r==1
    ylabel('Filter (dB)')
    end
    xlabel('f/fs')
end

for r=1:n   % set gray scales of last row of plots
    DG_SetTraceTint(30,2,FigNum,r)
    DG_SetTraceWidth(1,'*',FigNum,r)
end

for r=1:n   % set gray scales of last row of plots
    DG_SetTraceTint(30,2,FigNum,r+3)
    DG_SetTraceWidth(1,'*',FigNum,r+3)
end

% export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigIVnr2b.pdf', gcf);
        
