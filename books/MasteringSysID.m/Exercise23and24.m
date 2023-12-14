% Chapter 2 Exercise 23 and 24
% Generation of a multisine with a reduced crest factor
% Exc. 12: using random phase generations
% Exc. 13: using a crest factor minimization algorithm
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters

Ndata1=256;                      % period length
NTrial=10000;                    % number of trials for the random phase realizations
t=[0:Ndata1-1]';

% Part 1: random phase realizations
uMAll=zeros(NTrial,1);                      % initialization
ExcitationLines=1+floor([1:0.1*Ndata1]');   % excited lines of the multisine
Amp=ones(size(ExcitationLines)); % set the amplitudes
uMax=inf;                                   % initialization

for k=1:NTrial  % run over the Ntrial realizations    
    
	U=zeros(Ndata1,1);
	U(ExcitationLines)=Amp.*exp(j*2*pi*rand(length(ExcitationLines),1));
	u=2*real(ifft(U));u=u/std(u);           % generate 1 period and normalize
	if max(abs(u))<uMax                     % select the signal with smallest peak
        uSmall=u;
        uMax=max(abs(u));
	end
    
    uMaxAll(k)=uMax;
end

% Part 2: Crest factor minimization

UMinCrest=crestmin(fiddata([],Amp,ExcitationLines-1));      % FDiDENT toolbox call
z=msinprep(UMinCrest,Ndata1);                               % FDiDENT toolbox call
uMinCrest=z.input;uMinCrest=uMinCrest/std(uMinCrest);       % normalize

% Plot the results
FigNum=10                       % plot both results
figure(FigNum)
clf

subplot(2,2,1)
plot(t,uSmall,'k');
axis([0 Ndata1 -4 4])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)
title('Best random multisine')

subplot(2,2,2)
plot(t,uMinCrest,'k')
axis([0 Ndata1 -4 4])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)
title('Optimized multisine')

subplot(2,2,3)
edges=[-3:0.05:3];
H=histc(uSmall,edges);
bar(edges,H,'k')
axis([-4 4 0 15])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)
xlabel('u'),ylabel('Counts')

subplot(2,2,4)
H=histc(uMinCrest,edges);
bar(edges,H,'k')
axis([-4 4 0 15])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetLineWidth(1.0)
xlabel('u(k)'),ylabel('Counts')

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('MultisineSmallCrest1a.pdf', gcf);

FigNum=11                       % evolution of the peak value 
figure(FigNum)
semilogx([1:NTrial],uMaxAll','k',[1:NTrial],max(abs(uMinCrest))*ones(size(uMaxAll)),'--k')
axis([0 NTrial 0 4])
DG_SetFontSize(12)
DG_SetTraceWidth(0.5,1,FigNum)
DG_SetTraceWidth(1,2,FigNum)
DG_SetTraceGray(0.6,2,FigNum)
DG_SetLineWidth(1.0)
title('Evolution of the peak value')
xlabel('Trial number')

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('MultisineSmallCrest1b.pdf', gcf);


