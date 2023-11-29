% Chapter 2 Exercise 32
% Amplitude distribution of filtered noise
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all

% define the input parameters
Ndata=100000;                   % length data record
CutOff=0.1                      % Cut off frequency, relatively to fSample
[b,a]=butter(6,CutOff*2);       % noise shaping filter


u=sign(randn(Ndata,1));         % binaray input
y=filter(b,a,u);                % filtered output   

% analyse and plot the result
FigNum=5;
figure(FigNum)
clf

subplot(1,2,1)
edges=[-1.5:0.1:1.5];           % set HIST parameters
H=histc(u,edges);
bar(edges,H,'k')
title('Binary input')
xlabel('Input')
ylabel('Amplitude')
DG_SetFontSize(10)
DG_SetLineWidth(1.0)

subplot(1,2,2)
edges=[-1.5:0.1:1.5];
H=histc(y,edges);
bar(edges,H,'k')
title('Filtered binary input')
xlabel('Input')
ylabel('Amplitude')
DG_SetFontSize(10)
DG_SetLineWidth(1.0)

% Export the plot
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('Noise5.pdf', gcf); 