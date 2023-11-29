% Chapter 1 Exercise 13
% Noise on input and output: the errors-in-variables method 
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 25 November 2010

clear all
N=5000              % number of datapoints
NRep=100000         % number of repeated experiments
R0=1000             % exact value resistor
stdi0=0.01          % std. dev. input current
stdi=0.001          % std. dev. input current noise
stdu=1              % std. dev. voltage noise


i0=stdi0*randn(N,1);% generate the exact current

u0=R0*i0;           % exact voltage


vari=stdi^2;varu=stdu^2;            % variances input/output noise

for s=NRep:-1:1                     % loop over the different realizations
    if rem(s,1000)==0,s,end         % display the loop index
    
    iNoise=stdi*randn(N,1);         % white disturbing noise on current 
    uNoise=stdu*randn(N,1);         % white disturbing noise on voltage
    
    i=i0+iNoise;                    % add noise to current
    u=u0+uNoise;                    % add noise to voltage
    
    RLS(s,1)=(i'*u)/(i'*i);         % Least squares estimate
    
    
    z1=(u'*u/varu-i'*i/vari);       % calculate errors-in-variables estimate
    z2=u'*i/varu;
    REV(s,1)=(z1+sqrt(z1^2+4*z2^2/vari))/(2*z2);
    
end

BinWidth=0.5;X=[980:BinWidth:1010];     % pdf estimation parameters

pdfRLS=hist(RLS,X)/BinWidth/NRep;       % estimated pdf EIV
pdfIV=hist(REV,X)/BinWidth/NRep;        % estimated pdf IV

% plot the results
FigNum=5;
figure(FigNum),clf

plot(X,pdfRLS,'k',X,pdfIV,'k')
axis([980 1010 0 0.25])
    
DG_SetFontSize(12)
ylabel('pdf(R)')
xlabel('R')
DG_SetTraceTint(30,2,FigNum,1)
DG_SetTraceWidth(1,'*',FigNum,1)

% export the plot  
% DG_Init4PDF(gcf, 4.5);        % fixing the size, half standard height
% DG_MakePDF('FigEIVnr1.pdf', gcf); 
      


       
    


