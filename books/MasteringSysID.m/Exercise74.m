% Chapter 4 Exercise 74
% Identification in the frequency domain using nonparametric noise models
% and a random excitation
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% Define the input variables
N=5000;                             % length of the processed data records
NTrans=500;                         % eliminate transients of the simulation
NRep=100;                           % number of repeated simulations
IterMax=100;                        % maximum number of itereations
OrderG=2;                    	    % order of the plant
fMax=0.4;                           % to set the bandwidth of the systems (fs=1)
stdNoise=0.1;                       % standard deviation of disturbing noise
fMaxPlot=N*min(fMax*1.5,0.5);       % set frequencies for the plot
w=[1:fMaxPlot/100:fMaxPlot]/N*2*pi; %freq grid

% Generator filter
[bGen,aGen]=butter(3,2*fMax);       % filter used to generate the excitation of the system
GGen=freqz(bGen,aGen,w);            % FRF of the generator filter

% G0
[b0,a0]=cheby1(OrderG,10,2*fMax*0.9); % plant to be estimated
G0=freqz(b0,a0,w);                  % FRF of G0

% noise models
OrderNoise=2;                       % order of the noise filter
[bNoise,aNoise]=butter(OrderNoise,2*fMax/3); % noise model
GNoise=freqz(bNoise,aNoise,w);      % FRF noise filter


% Generate and process the data
counter=0;   % count the number of converged runs

for s=NRep:-1:1      % loop over de realizations
      
    u0=randn(N+NTrans,1);u0=filter(bGen,aGen,u0);       % filtered input noise
    y0=filter(b0,a0,u0);                                % exact output
    yNoise=filter(bNoise,aNoise,randn(N+NTrans,1));yNoise=stdNoise*yNoise;
                                                        % disturbing noise
    y=y0+yNoise;                                        % disturbed output
    u0(1:NTrans)=[];y(1:NTrans)=[];                     % eliminate transients of the simulation
    
    % estimate Box-Jenkins model
    Mbj= bj([y u0],[OrderG+1 OrderNoise OrderNoise OrderG 0],'maxiter',IterMax,'tol',1e-6,'lim',0);
 
    % extract results if converged
    if max([Mbj.EstimationInfo.Iterations ])<IterMax    % converged?
        counter=counter+1;
        [s counter]                                     % follow the simulation                           
        % Box Jenkins processing
        [b,a]=tfdata(Mbj,'v');            
        G=freqz(b,a,w);
        Gbj(:,counter)=G;
        ebj(:,counter)=abs(G0-G);
     end  % end if loop
    
%        
% Non Parametric noise model with FD-identification
%

% processing with the local polynomial method
    U0=fft(u0(:))/sqrt(N);U0(N/2+1:end)=[];                 % FFT input
    Y=fft(y(:))/sqrt(N);Y(N/2+1:end)=[];                    % FFT output

    Data.U = U0(:).';                                       % store spectra as rows
    Data.Y = Y(:).';                                        % store spectra as rows
    Data.Freq=[0:N/2]'/N;                                   % store corresponding frequencies

    Method.moment=6;                                        % To have enough degrees of freedom in the result
    [CY, Ym, TY, GLocPol1, CvecG2OE,M] = LocalPolyAnalv3(Data,Method);
                                                            % Estimate nonparametric model
    YLocPolVar=squeeze(CY.m_nt);                            % the estimated variance
    YLocPol=squeeze(Ym.m_nt);                               % the estimated output after transient removal 

    fs=1;fLocPol=[0:N/2-1]'/N*fs;                           % create freq. vector Loc. Pol. Method

    % Frequency Domain identification
    Fdat=fiddata(YLocPol(:),U0,fLocPol,YLocPolVar,0);       % create Frequency Domain Data object
    nb=OrderG;      % order numerator
    na=OrderG;      % order denominator
    MFD=elis(Fdat,'z',nb,na,struct('fs',1));                % frequency domain identification
    [domain,bFD,aFD]=imppar(MFD);                           % processing the estimates
    GFD(:,s)=freqz(bFD,aFD,w);
   
end  % end loop over the repetitions

f=w/2/pi;
    
%Calculate mean square error
for k=length(f):-1:1
RMSbj(k,1)=sqrt(mean(abs(Gbj(k,:)-G0(k)).^2));
RMSGFD(k,1)=sqrt(mean(abs(GFD(k,:)-G0(k)).^2));
end
stdbj=std(Gbj,0,2);

% Plot the results
FigNum=1;
figure(FigNum)
clf

plot(f,db(G0),'k',f,db(RMSbj),'k',f,db(RMSGFD),'k')
axis([0 fMax -80 10])
ylabel('Amplitude (dB)')
xlabel('Frequency');
DG_SetFontSize(18,FigNum)        
DG_SetTraceTint(100,1,FigNum)     
DG_SetTraceWidth(2,1,FigNum)   
DG_SetTraceTint(100,2,FigNum)     
DG_SetTraceWidth(0.5,2,FigNum)    
DG_SetTraceTint(50,3,FigNum)    

% Export the plot
%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('FigNonParNoiseModel.pdf', gcf);   