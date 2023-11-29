% Chapter 7 Exercise 98
% Variance of the parametric estimate of the BLA of y=u0^n
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 8 December 2010


% Exercise 98 part 1: a cubing system y=u.^3        varying the number of data points
% Exercise 98 part 2: a nth degree system y=u.^n    varying degree n

clear all

%
% Exercise 98 part 1:  y=u.^3
%

NAll=[10 100 1000 10000]        % number of points per test
NRep=100                        % number of repeated tests
n=3                             % degree of the nonlinear system
x=[1:NRep]';                    % x-axis for plots

for r=length(NAll):-1:1         % run over the different record lenghts
    N=NAll(r)                   % select the number of points
    for k=NRep:-1:1                 % run over the repeated experiments 
        u=randn(N,1);               % white gaussian input
        y=u.^n;                     % 3th degree NL    
        a(k,r)=(y'*u)/(u'*u);       % LS-estimate BLA
        e=y-a(k,r)*u;               % error 
        stdTh(k,r)=sqrt(var(e)/mean(u.^2)/N);  % theoretic std. dev. estim in Linear framework   
    end
end    % end loop over number of points per test 

FigNum=1
figure(FigNum),clf

for r=1:4
subplot(2,2,r)
plot(x,std(a(:,r))*sqrt(NAll(r)),'k',x,stdTh(:,r)*sqrt(NAll(r)),'k')   % plot the normalized std. dev.
axis([0 100 0 8])
title(['N = ',num2str(NAll(r))]) 
if r>2,xlabel('Exp. Nr.'), end
if or(r==1,r==3), ylabel('Std. Dev.'),end
end

DG_SetFontSize(18,FigNum,'*')         
DG_SetTraceWidth(2,'*',FigNum,'*')  %


%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('StdaBLAUnderEstimate.pdf', gcf);   

%
% Exercise 98 Part 2:  y=u.^n
%

clear all
nAll=[  3  5  7  9];            % degree of the NL
N=1000;                         %number of points per test
NRep=10000                      % number of repeated tests to estimate the experimental variance
x=[1:NRep]';                    % x-axis for plots

for r=length(nAll):-1:1         % run over the different record lenghts
    n=nAll(r)   
    for k=NRep:-1:1             % run over the repeated experiments
        u=randn(N,1);           % white gaussian input
        y=u.^n;                 % 3th degree NL
        a(k)=(y'*u)/(u'*u);     % BLA
        e=y-a(k)*u;             % error 
        varTh(k)=(var(e)/var(u)/N);  % theoretic std. dev. estim in Linear framework
    end                         % end loop over the realizations
    stdExp(r)=std(a);
    stdTh(r)=sqrt(mean(varTh));
end                             % end loop over degrees of NL 

FigNum=2
figure(FigNum),clf

subplot(1,2,1)
plot(nAll,db(stdExp),'+k') 
axis([0 10 -30 50])
xlabel('Degree NL')
ylabel('Std. Dev. (dB)')

subplot(1,2,2)
plot(nAll,db(stdExp./stdTh),'+k',nAll,db(2*nAll+1)/2,'--k')
axis([0 10 0 14])
xlabel('Degree NL')
ylabel('Error Std. Dev. (dB)')

DG_SetFontSize(18,FigNum,'*')         
DG_SetTraceWidth(2,'*',FigNum,'*')  % 

  

%DG_Init4PDF(gcf, 6);        % fixing the size, half standard height
%DG_MakePDF('StdaBLAxn.pdf', gcf);   

