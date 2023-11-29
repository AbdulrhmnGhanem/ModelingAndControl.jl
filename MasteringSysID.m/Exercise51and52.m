function Excercise51and52

% Chapter 3 Exercise 51 and Exercise 52
% Exercise 51: The local polynomial method
% Exercise 52: Estimation the power spectrum of the disturbing noise
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 29 November 2010


clear all

% define the input variables
M=16;                   % number of blocks to be processed
NPeriod=128             % length of a block
N=NPeriod*M;            % Number of data points in the measurement
fs=NPeriod              % sampling frequency
Ampl=1                  % RMS-value of the multi sine
NTrans=1024             % eliminate transient effects simuation

t=[0:NPeriod-1]/fs;     % time vector 1 period
tAll=[0:N-1]/fs;        % time vector full sequence
Lines=[1:NPeriod/2];    % select the lines to be processed and plotted
f=[Lines-1]'/NPeriod*fs; % frequency vector

[b,a]=cheby1(2,10,0.2);             % test system

% generate excitation signal
u=randn(N+NTrans,1)*Ampl;
    
% system response    
y=filter(b,a,u);                    % response of the system  
y(1:NTrans)=[];u(1:NTrans)=[];      % eliminate transients
    
% processing using the hanning window    
y=reshape(y,NPeriod,M);             % one block per collumn
u=reshape(u,NPeriod,M);
wWind=hann(NPeriod,'periodic');wWind=kron(wWind,ones(1,M)); % prepare the window function
Y=fft(y.*wWind)/sqrt(NPeriod/2);Y=Y(Lines,:);               % spectral analysis, block per block
U=fft(u.*wWind)/sqrt(NPeriod/2);U=U(Lines,:);

UU=mean(abs(U(:,1:M)).^2,2);        % SUU with Hanning window
YU=mean(Y(:,1:M).*conj(U(:,1:M)),2);% SYU with Hanning window
G=YU(:)./UU(:);                     % FRF estimate with Hanning window

G0Hann=freqz(b,a,f/fs*2*pi);        % true FRF of the system at Hann frequencies
    
% processing with the local polynomial method: prepare the data + call
U=fft(u(:))/sqrt(N);Y=fft(y(:))/sqrt(N);U(N/2+1:end)=[];Y(N/2+1:end)=[];
                                    % FFT analysis
Data.U = U(:).';                    % store the spectra as raws
Data.Y = Y(:).';                    % store spectra as rows
Data.Freq=[0:N/2]'/N;               % store corresponding frequencies

method.order=2;                     % 2nd degree polynomial
method.bandwidth=6;                 % width of the window
method.transient=1;                 % =1 --> estimate the transients

[GLocPol varY]=LocalPolynMethod(Data,method);
                                    % apply local polynomial method
fLocPol=[0:N/2-1]'/N*fs;            % create freq. vector Loc. Pol. Method
G0LocPol=freqz(b,a,fLocPol/fs*2*pi);% true FRF of the system at Hann frequencies  

% plot the results
FigNum=1;
figure(FigNum) % plot estimated FRF
clf
plot(f,db(G),'k',f,db(G-G0Hann),'k.',fLocPol,db(GLocPol-G0LocPol),'k.')
axis([0 65 -150 5])
ylabel('Amplitude (dB)')
xlabel('f (Hz)')

DG_SetTraceTint(30,3,FigNum,'*')
DG_SetTraceWidth(1,'*',FigNum,'*')

% export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigExcitationLocPol1.pdf', gcf);

FigNum=2
figure(2)  % plot noise estimate
clf
plot(fLocPol,db(Y),'.k',fLocPol,db(varY)/2,'.k')

axis([0 65 -150 5])
ylabel('Amplitude (dB)')
xlabel('f (Hz)')

DG_SetTraceTint(30,2,FigNum,'*')
DG_SetTraceWidth(1,'*',FigNum,'*')

% export the plot
% DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
% DG_MakePDF('FigExcitationLocPol2.pdf', gcf);



 function [G,varY]=LocalPolynMethod(DataPoly,method)

% Local polynomial method
% f: frequency vector
% G: estimated FRF at f

% DataPoly.U    input Fourier coefficients
%         .Y    output Fourier coefficients
%
% method.bandwidth   the bandwidth in which the solution is calculated
% method.order       the degree of the polynomial
% method.transient   determines the estimation of the transient term (optional; default 1)  
%                                                   1: transient term is estimates 
%                                                   0: no transient term is estimated 


w=method.bandwidth;
U=DataPoly.U(:);                    % input DFT
Y=DataPoly.Y(:);                    % output DFT
N=length(U);

% window lines
w1=ceil(w/2);                       % calculation of the window width
wLines=[-w1:w1];

for k=N:-1:1 % frequency loop
    SelectLines=k+wLines;           % select the lines to be used at that frequency
    kFreq=k;                        % frequency at which an estimate is made
    fFit=wLines;
    % edges setting
    if SelectLines(end)>N;          % keep window within upperlimits of the measurement band
        fShift=SelectLines(end)-N;
        SelectLines=SelectLines-fShift;
        fFit=wLines-fShift;
    end
                     
    if SelectLines(1)<1;            % keep window within lowerlimits of the measurement band
        fShift=abs(SelectLines(1))+1;
        SelectLines=SelectLines+fShift;
        fFit=wLines+fShift;
    end   % end edges setting
    
    % the polynomial fit
    UB=U(SelectLines);               % select lines to be fitted
    YB=Y(SelectLines);
    fFit=SelectLines-k;              % select the corresponding freq
    
    fFitScaled=fFit/w1;              % normalize the frequencies between [-1 1]

    [Gfit,Yvar]=LocalLPMLocal(YB,UB,fFitScaled,method);
    
    % evaluate the results
    G(k,1)=Gfit(1);
    varY(k,1)=Yvar;
     
end   % end frequency loop



function [G,Yvar]=LocalLPMLocal(Y,U,f,method)
Y=Y(:);
U=U(:);
f=f(:);                  % complex frequency

nOrder=method.order;     % polynomial in f

if method.transient==1   % transient estimation requested?
    select=1;
else
    select=0;
end                 % estimate transients

K=zeros(length(f),(nOrder+1)*(1+select));

K(:,1)=U;counter=1;
for k=1:nOrder      % U.*f^k    input contributions
    counter=counter+1;
    K(:,counter)=U.*(f.^k);
end

if select==1
    for k=0:nOrder  % 1.*f^k    transient contributions
        counter=counter+1;
        K(:,counter)=f.^k;    
    end
end                  % end transient contribution

% solve equations
Theta=K\Y;

G=Theta(1);
p=Theta(1:nOrder+1).';


% estimate the variance
Yerror=Y-K*Theta;              % calculate residual
q=(length(Y)-length(Theta));   % number of (complex) degrees of freedom
Yvar=sum(abs(Yerror).^2)/q;    % estimate complex variance 

