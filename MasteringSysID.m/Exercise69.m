% Chapter 4 Exercise 69
% Generating uncertainty bounds for estimated models
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 2 December 2010

clear all

% define the input variables
N=1000;                                 % length of the processed data records
NTrans=500;                             % eliminate transients of the simulation
NRep=1000;                              % number of repeated simulations
IterMax=100;                            % maximum number of itereations

OrderG=2;                               % order of the plant
fMax=0.2;                               % max frequency of the generation filter for random excitation
[bGen,aGen]=butter(3,2*fMax*1.5);       % filter used to generate the excitation of the system
stdNoise=0.1;                           % standard deviation of disturbing noise

% G0
b0=[0.2 0.6 0.3];                       %a direct term is present! b0 not equal to zero!
a0=[1 -1.2 0.9];
w=[1:N/2]'/N*2*pi;                      %frequency grid for the plot
G0=freqz(b0,a0,w);                      % FRF exact system

% noise model
OrderNoise=2;            
bNoise=[1  0.7  0.7]; 
aNoise=[1 -0.8  0.3];
GNoise=freqz(bNoise,aNoise,w);          % FRF of the noise model

counter=0;                              % count the number of converged runs

for s=NRep:-1:1                         % loop over de realizations
      
    u0=randn(N+NTrans,1);u0=filter(bGen,aGen,u0);   % filtered random input
    y0=filter(b0,a0,u0);                            % exact output
    yNoise=filter(bNoise,aNoise,randn(N+NTrans,1)); % disturbing noise
    yNoise=stdNoise*yNoise;                         % scaling the noise
    y=y0+yNoise;                                    % disturbed output
    u0(1:NTrans)=[];y(1:NTrans)=[];                 % eliminate transients of the simulation

    Mbj{s}= bj([y u0],[OrderG+1 OrderNoise OrderNoise OrderG 0],'maxiter',IterMax,'tol',1e-6,'lim',0);
                                                    % Box-Jenkins estimation
     
    % extract results if converged
 
    if ([Mbj{s}.EstimationInfo.Iterations])<IterMax   % converged?
                           counter=counter+1;
        [s counter]                       
        % Box Jenkins processing
        [b,a]=tfdata(Mbj{s},'v');            
        G=freqz(b,a,w);
        Gbj(:,counter)=G;
        ebj(:,counter)=abs(G0-G); 
    end  % end if loop
    
end  % end loop over the repetitions
    
%Calculate standard deviation
STDbj=std(abs(Gbj),0,2);                            % from simulations
STDbjPh=std(angle(Gbj),0,2);                        % from simulations

   
% plot the results
FigNum=17
figure(FigNum),
clf

subplot(1,2,1)                                      % amplitude results
s=1;
[mag,phse,wbj,STDbjTh,sdphase]=bode(Mbj{s},w);STDbjTh=STDbjTh(:);
                                % from model

f=w/2/pi;
plot(f,db(G0),'k',f,db(STDbj),'k',f,db(STDbjTh),'k')
axis([0 0.5 -55 25])
xlabel('Frequency');
ylabel('Amplitude');
DG_SetFontSize(12,FigNum,1)
DG_SetTraceTint(100,1,FigNum,1)    
DG_SetTraceWidth(1,1,FigNum,1)   
DG_SetTraceTint(30,2,FigNum,1)    
DG_SetTraceWidth(3,2,FigNum,1)    
DG_SetTraceTint(100,3,FigNum,1)    
DG_SetTraceWidth(1,3,FigNum,1)     



subplot(1,2,2)                                      % phase results 
[mag,phse,wbj,STDbjTh,sdphaseTh]=bode(Mbj{s},w);sdphaseTh=sdphaseTh(:);
                                % from model

plot(f,STDbjPh*180/pi,'k',f,sdphaseTh,'k')
axis([0 0.5 0 1.5])
xlabel('Frequency');
ylabel('Phase (deg)');
DG_SetFontSize(7,FigNum,2)
DG_SetTraceTint(30,1,FigNum,2)    % simulated std. dev.
DG_SetTraceWidth(1.5,1,FigNum,2)   % 
DG_SetTraceTint(100,2,FigNum,2)    % theoretic std. dev. 
DG_SetTraceWidth(0.5,2,FigNum,2)     % 

% Export the results
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigBJStdDev.pdf', gcf);  

