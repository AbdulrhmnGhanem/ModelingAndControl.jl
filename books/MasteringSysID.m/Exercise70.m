% Chapter 4 Exercise 70
% Study of the behavior of the BJ-model in combination with prefiltering
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
NRep=100                    % number of simulations
N=1024;                     % data length
NTrans=1024;                % eliminate the transient of the simulation
stdNoise=0.1;               % standard deviation of disturbing noise
OrderG=2;                   % order of the plant
OrderNoise=2;               % BJ noise model
fGen=0.25;                  % bandwidth of the generator signal (fs=1)
IterMax=100;                % maximum number of iterations in BJ

[bGen,aGen]=butter(3,2*fGen);bGen(2)=0.9*bGen(2);       % filter to generate to random excitation

[b0,a0]=cheby1(OrderG,5,2*0.08);b0(2)=b0(2)*1.3;        % plant model

[b,a]=butter(OrderNoise,2*0.2); b=b+0.1*a;              % noise model
bNoise=b;   
aNoise=a; 


fMaxPlot=N*0.5;                                         % set up frequency grid for plots
w=[0:fMaxPlot/100:fMaxPlot]/N*2*pi; f1=w/2/pi;          % freq grid for plot
G0=freqz(b0,a0,w);                                      % FRF plant

OrderPF=3;                                              % define the prefilter
fmax=0.1;                                               % frequency band of the prefilter
[bf,af]=butter(OrderPF,fmax*2);bf=bf+0.01*af;           % create lowpss prefilter
OrderPFBJ=OrderPF;                                      % additional order for BJ-filt noise model
GPF=freqz(bf,af,w);                                     % FRF prefilter

% generate and process the data
counter=0
for s=NRep:-1:1                                         % loop over de realizations
    [s]                                                 % follow the simulation
      
    % generate the signals, no prefiltering
    u0=filter(bGen,aGen,randn(N+NTrans,1));             % generate colored noise excitation
    y0=filter(b0,a0,u0);                                % exact output
    yNoise=filter(bNoise,aNoise,randn(N+NTrans,1));yNoise=stdNoise*yNoise;
                                                        % disturbing noise
    y=y0+yNoise;                                        % disturbed output
   
    u0PF=filter(bf,af,u0);                              % prefiltered input
    yPF=filter(bf,af,y);                                % prefiltered output

    % process, no prefiltering
    u0(1:NTrans)=[];y(1:NTrans)=[];                     % eliminate transients simulation
    u0=u0(:);y=y(:);                                    % make colum vectors
  
    % prefiltered signals
    u0PF(1:NTrans)=[];yPF(1:NTrans)=[];                 % eliminate transients simulation
    u0PF=u0PF(:);yPF=yPF(:);                            % make colum vectors

% Box-Jenkins estimation
    MBJ= bj([y u0],[OrderG+1 OrderNoise OrderNoise OrderG 0],'maxiter',1000,'tol',1e-4,'lim',0);
                                                        % no prefiltering
    MBJPF= bj([yPF u0PF],[OrderG+1 OrderNoise OrderNoise OrderG 0],'maxiter',1000,'tol',1e-4,'lim',0);
                                                        % prefiltering, no increased noise model
    MBJPF2= bj([yPF u0PF],[OrderG+1 OrderNoise+OrderPFBJ OrderNoise+OrderPFBJ OrderG 0],'maxiter',1000,'tol',1e-4,'lim',0);
                                                        % prefiltering, increase noise model

% extract results
    
     if max([MBJ.EstimationInfo.Iterations MBJPF.EstimationInfo.Iterations ,...
                               MBJPF2.EstimationInfo.Iterations])<IterMax   % converged?
        counter=counter+1

        [bBJ,aBJ]=tfdata(MBJ,'v');            
        G=freqz(bBJ,aBJ,w);
        GBJ(:,counter)=G;
        eBJ(:,counter)=(G0-G);

        [bBJPF,aBJPF]=tfdata(MBJPF,'v');            
        G=freqz(bBJPF,aBJPF,w);
        GBJPF(:,counter)=G;
        eBJPF(:,counter)=(G0-G);

        [bBJPF2,aBJPF2]=tfdata(MBJPF2,'v');            
        G=freqz(bBJPF2,aBJPF2,w);
        GBJPF2(:,counter)=G;
        eBJPF2(:,counter)=(G0-G);
   
     end
    
end  % end loop over the repetitions
    
%Calculate mean square error
EBJ=sqrt(mean(abs(eBJ).^2,2));
EBJPF=sqrt(mean(abs(eBJPF).^2,2));
EBJPF2=sqrt(mean(abs(eBJPF2).^2,2));

% calculate the std. dev.
Gstd=std(abs(GBJ),0,2);                 % std. on the magnitude
GstdPF=std(abs(GBJPF),0,2);             % std. on the magnitude
GstdPF2=std(abs(GBJPF2),0,2);           % std. on the magnitude


% plot the results 
FigNum=18;
figure(FigNum)
clf
plot(f1,db(G0),'k',f1,db(GPF),'k',f1,db(EBJ),'k',f1,db(EBJPF),'k',f1,db(EBJPF2),'k')
axis([0 0.5 -60 5])
ylabel('Amplitude (dB)');
xlabel('Frequency');
DG_SetFontSize(7,FigNum)
DG_SetTraceTint(100,1,FigNum)     
DG_SetTraceWidth(1.25,1,FigNum)     
DG_SetTraceTint(100,2,FigNum)     
DG_SetTraceWidth(0.5,2,FigNum)       
DG_SetTraceTint(100,3,FigNum)     
DG_SetTraceWidth(0.5,3,FigNum)      
DG_SetTraceTint(50,4,FigNum)     
DG_SetTraceWidth(1,4,FigNum)       
DG_SetTraceTint(100,5,FigNum)     
DG_SetTraceWidth(0.5,5,FigNum)       
DG_SetTraceStyle(':',5,FigNum)

% Export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigBJPrefilter1.pdf', gcf);  


%
% verify the actual observed oncertainties with the BJ-uncertainty bounds
% use the model of the last simulation
%

% experimental std. dev. on the amplitude
FigNum=19;
figure(FigNum)
clf

w=f1*2*pi;
[mag,phse,wbj,sdmag,sdphase]=bode(MBJ,w);mag=mag(:);sdmag=sdmag(:);
[mag,phse,wbj,sdmagPF,sdphase]=bode(MBJPF,w);mag=mag(:);sdmagPF=sdmagPF(:);
[mag,phse,wbj,sdmagPF2,sdphase]=bode(MBJPF2,w);mag=mag(:);sdmagPF2=sdmagPF2(:);

plot(f1,db(G0),'k',f1,db(sdmagPF2),'k',f1,db(GstdPF2),'k--',...
                       f1,db(sdmagPF),'k',f1,db(GstdPF),'k--'),shg
axis([0 0.5 -60 5])
ylabel('Amplitude (dB)');
xlabel('Frequency');
DG_SetFontSize(7,FigNum)
DG_SetTraceTint(100,1,FigNum)    
DG_SetTraceWidth(1.25,1,FigNum)   
DG_SetTraceTint(100,2,FigNum)    
DG_SetTraceWidth(0.5,2,FigNum)      
DG_SetTraceTint(100,3,FigNum)    
DG_SetTraceWidth(0.5,3,FigNum)    
DG_SetTraceTint(50,4,FigNum)    
DG_SetTraceWidth(0.5,4,FigNum)    
DG_SetTraceTint(50,5,FigNum)    
DG_SetTraceWidth(0.5,5,FigNum)     

% Export the plot
% DG_Init4PDF(gcf, 6);        % 
% DG_MakePDF('FigBJPrefilter2.pdf', gcf);  