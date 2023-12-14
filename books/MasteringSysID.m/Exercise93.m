function Exercise93

% Chapter 6 Exercise 93
% Prediction of the bias contribution in the BLA
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 9 December 2010


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %% 
% Prediction of the difference between the BLA and the true   %%
% underlying linear system via the odd stochastic nonlinear   %%
% distortions. No noise is added to the input/output signals. %%
%                                                             %%
% The system:                                                 %% 
% Wiener-Hammerstein system, consisting of the cascade of a   %%
% second order discrete-time system, a polynomial static      %%
% nonlinearity, and a second order discrete-time system.      %%
%                                                             %%
% The signal:                                                 %% 
% Odd multisine with random harmonic grid.                    %%
%                                                             %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition excitation signal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 700e6;
fmax = 1100e6;
fres = 1e6;
fs = 4e9;
frat = 1.02;
Nblock = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice frequency spacing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spacing = 'lin';
% Spacing = 'log'; fmin = 10e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice type multisine  %
%%%%%%%%%%%%%%%%%%%%%%%%%%

TypeMultisine = 'odd';
urms = 1;       % rms value multsine excitation


DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fs = fs;
DefFreq.fres = fres;
DefFreq.frat = frat;
[ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMultisine);

% odd multisine with all odd harmonics excited
StarHar = round(fmin/(fs/N));
StopHar = round(fmax/(fs/N));
ExcitedHarm = ([StarHar+1:2:StopHar-1]).';


%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation measurement %
%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 400;                     % number of realisations random phase multisine
% M = 1;
P = 2;                      % number of consecutive periods

freqall = ([0:1:N-1]*fs/N).';               % all DFT frequencies over full unit circle
qall = exp(-sqrt(-1)*2*pi*freqall/fs);      % z^-1 over full unit circle
freq = ExcitedHarm*fs/N;                    % the excited frequencies
qOdE = exp(-sqrt(-1)*2*pi*freq/fs);   % z^-1 at odd excited frequencies
fstart = 600e6;                             % first measured frequency
fstart = 1e6;
if Spacing == 'log'
    if strcmp(TypeMultisine, 'odd')
        fstart = fres / 2;
    else % full multisine
        fstart = fres ;
    end % if
    fstop = fmax * 3;
end % if
fstop = 1200e6;                             % last measured frequency
fstop = 2000e6;                         % last measured frequency
if strcmp(TypeMultisine, 'odd')
    freqmeas = ([fstart:fres/2:fstop]).';       % all measured frequencies
else % full multisine
    freqmeas = ([fstart:fres:fstop]).';       % all measured frequencies
end % if
F = length(ExcitedHarm);                    % number of measured frequencies
Uall = zeros(M, P, F);                   % all input spectra at the measured frequencies
Yall = zeros(M, P, F);                   % all output spectra at the measured frequencies


% discrete-time LTI1
a1 = [1, -0.2, 0.9];
b1 = [1];
G1all = polyval(fliplr(b1), qall)./polyval(fliplr(a1), qall);
G1 = polyval(fliplr(b1), qOdE)./polyval(fliplr(a1), qOdE);

% polynomial static nonlinearity
alfa = 0.1;                                % coefficient x^2
beta = 0.001;                              % coefficient x^3
if Spacing == 'log'
    beta = 0.1;
end % if

% discrete-time LTI2
a2 = [1, -0.5, 0.9];
b2 = [1, 0.5];
G2all = polyval(fliplr(b2), qall)./polyval(fliplr(a2), qall);
G2 = polyval(fliplr(b2), qOdE)./polyval(fliplr(a2), qOdE);


for ii=1:M
    
    % generation one period random phase multisine with rms value one
    u = CalcMultisine(ExcitedHarm, N);
    u = u*urms;
    
   % response to the first LTI system
    U = fft(u)/sqrt(N);
    X = U.*G1all;
    x = real(ifft(X)*sqrt(N));

    % response of the static nonlinearity
    z = x + alfa*x.^2 + beta*x.^3;

    % response to the second LTI system
    Z = fft(z)/sqrt(N);
    Y = Z.*G2all;
    
    % collect measurements
    Uall(ii,:,:) = repmat((U(ExcitedHarm+1)).', P, 1);
    Yall(ii,:,:) = repmat((Y(ExcitedHarm+1)).', P, 1);
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the mean FRF, its sample variance, the noise variance, % 
% and the variance of the stochastic nonlinear contributions       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[G] = BLA_per(Yall./Uall);

% true underlying linear system
G0 = G1.*G2;
% best linear approximation
GBLA = G.mean;
% stochastic nonlinear distortions
Gs = G.stds;
% bias contribution
Gb = G0-GBLA;

figure(1);
freq = freq/1e9;    % frequency in GHz
plot(freq, db(GBLA), 'k', freq, db(Gs), 'b', freq, db(Gb), 'r')
xlabel('Frequency (GHz)')
ylabel('Amplitude (dB)')
title('G_{BLA}: black; G_S: blue; G_B: red')

function [G, Y, U, CYU] = BLA_per(Yall, Uall, Rall);
%
%		function [G, Y, U, CYU] = BLA_per(Yall, Uall, Rall);
%
%       
%           Nonlinear analysis of FRF or input/output measurements at F frequencies
%           made with different M phase realisations of random phase multisines with the
%           same rms value and amplitude spectrum. P consecutive periods are measured
%
%           Yall(m,p,k) = Y0[m,p](k) + NY[m, p](k) + YS[m](k)
%           Uall(m,p,k) = U0[m,p](k) + NU[m, p](k)
%           Rall(m,p,k) = R[m](k) = |R(k)| exp(sqrt(-1)*angle(R[m](k)))
%
%       Examples
%           1. [G] = BLA_per(Gall)                          where Gall is an M x P x F matrix of FRF's
%           2. [G] = BLA_perYall, Uall)                     where Yall, Uall are  M x P x F matrices of output, input spectra
%           3. [G, Y, U, CYU] = BLA_per(Yall, Uall, Rall)	where Yall, Uall are  M x P x F matrices of output, input spectra
%                                                         	and Rall is M x P x F or M x F matrix of the spectra of the reference signals
%
%	OUTPUT PARAMETERS
%
%	    G		:	struct('mean', 'stdNL', 'stdn', 'stds')
%                   G.mean  :   mean value of FRF estimate over the M x P FRF measurements
%                   G.stdNL :   total standard deviation mean FRF (contribution noise + NL distortion)
%                   G.stdn  :   noise standard deviation mean FRF
%                   G.stds  :   standard deviation stochastic NL distortions w.r.t. one multisine experiment
%                               = order of magnitude of NL bias errors = the difference between the best linear approximation G
%                                 and the true underlying linear system
%
%       Y       :   struct('mean', 'abs', 'stdn', 'stdNL', 'stds')
%                   Y.mean	:   if Rall is available then
%                                    Y.mean is the mean value of the projected output Yall./Rall
%                                    over all periods and realisations
%                               else
%                                   Y.mean is empty except if M = 1 then Y.mean is the mean value of the
%                                   output over all periods
%                   Y.abs	:   if Rall is available then
%                                   Y.abs = abs(Y.mean) 
%                               else
%                                   rms value |Y(k)| over different realisations compensated for the noise variance 
%                   Y.stdNL :   total variance Y.mean (as complex number) 
%                   Y.stdn  :   noise variance Y.mean (as complex number)
%                   Y.stds  :   output stochastic nonlinear distortions w.r.t. one multisine realisation 
%
%       U       :   struct('mean', 'abs', 'stdn', 'stdNL', 'stds')
%                   U.mean	:   if Rall is available then
%                                    U.mean is the mean value of the projected output Uall./Rall
%                                    over all periods and realisations
%                               else
%                                   U.mean is empty except if M = 1 then U.mean is the mean value of the
%                                   output over all periods
%                   U.abs	:   if Rall is available then
%                                   U.abs = abs(U.mean) 
%                               else
%                                   rms value |U(k)| over different realisations compensated for the noise variance 
%                   U.stdNL :   total variance U.mean (as complex number) 
%                   U.stdn  :   noise variance U.mean (as complex number)
%                   U.stds  :   input stochastic nonlinear distortions w.r.t. one multisine realisation 
%
%       CYU     :   struct('NL', 'n')
%                   CYU.NL  :   total covariance matrix of [Y; U], size 2 x 2 x F (noise + NL distortion)
%                                   CYU.NL(1, 1) = varY = E{|NY|^2}/(M*P) + E{|YS|^2}/M 
%                                   CYU.NL(2, 2) = varU = E{|NU|^2}/(M*P) + E{|US|^2}/M 
%                                   CYU.NL(1, 2) = varYU = E{NY * conj(NU)}/(M*P) 
%                                   CYU.NL(2, 1) = varUY = E{NU * conj(NY)}/(M*P) + E{YS * conj(US)}/M 
%                               with NX the noise on X w.r.t one period, and YS the stochastic NL distortions w.r.t. one experiment 
%                               Notes: 
%                                       - calculated only if Yall, Uall, and Rall are available 
%                                       - if the actuator is linear then US = 0 (no nonlinear interaction
%                                         between the generator and the NL system) 
%
%                   CYU.n   :   noise covariance matrix of [Y.mean; U.mean] as complex numbers, size 2 x 2 x F
%                                   CYU.n(1, 1) = varY = E{|NY|^2}/(M*P)
%                                   CYU.n(2, 2) = varU = E{|NU|^2}/(M*P)
%                                   CYU.n(1, 2) = varYU = E{NY * conj(NU)}/(M*P)
%                                   CYU.n(2, 1) = varUY = E{NU * conj(NY)}/(M*P) 
%                               with NX the noise on X w.r.t. one period
%
%	INPUT PARAMETERS
%
%	    Yall	:	M x P x F FRF or output spectrum measurements, where F is the number of frequencies,
%				    P the number of consecutive periods of one realisation of the random multisine excitation,
%                   and M is the number of independent realisations of the random phase multisine with
%                       1. the same rms value
%                       2. the same amplitude spectrum
%                   If only one argument is passed then Yall = Gall = FRF measurements
%       Uall    :   Optional argument, M x P x F input spectrum measurements
%       Rall    :   Optional argument, M x P x F or M x F matrix of the spectra of the M reference signals (signal stored in the waveform generator)
%
% Rik Pintelon, November 2004
% version December 2007
%


%%%%%%%%%%%%%%%%%%
% initialisation %
%%%%%%%%%%%%%%%%%%

[M, P, F] = size(Yall);
G = struct('mean', [], 'stdn', [], 'stdNL', [], 'stds', []);
Y = struct('mean', [], 'abs', [], 'stdn', [], 'stdNL', [], 'stds', []);
U = struct('mean', [], 'abs', [], 'stdn', [], 'stdNL', [], 'stds', []);
CYU = struct('NL', [], 'n', []);


%%%%%%%%%%%%%%%%%%%%%%%%%
% FRF measurements only %
%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 1)
    
	G.mean = squeeze(mean(mean(Yall, 2), 1));
	if M > 1
		G.stdNL = squeeze(std(mean(Yall, 2), 0, 1))/sqrt(M);
	end;
	if P > 1
		G.stdn = (squeeze(mean(std(Yall, 0, 2).^2, 1))/(M*P)).^0.5;
	end
    
end % if nargin = 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input/output measurements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2) | (nargin == 3)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If the reference signal is available then project      %
    % input and output on reference signal.                  %
    % It allows to average the spectra over the realisations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 3
       
        % calculating exp(j*angle(R))
        Rall = Rall ./ abs(Rall);
       
        % M x F matrix of reference spectra
        if length(size(Rall)) == 2
            dummy = zeros(M, 1, F);
            dummy(:, 1, :) = Rall;
            Rall = repmat(dummy, [1, P, 1]);
        end % if dim(Rall) = 2
        
        % turn input/output spectra with the phase of the reference signal
        Uall = Uall ./ Rall;
        Yall = Yall ./ Rall;
        
    end % if nargin = 3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mean and variance over P periods %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % averaging spectra over periods before FRF calculation
    Uallm = mean(Uall, 2);                                              % mean over consecutive periods
	Yallm = mean(Yall, 2);                                              % mean over consecutive periods
    if P > 1
		resUall = Uall - repmat(Uallm, [1, P, 1]);
		resYall = Yall - repmat(Yallm, [1, P, 1]);
		varUallm = squeeze((sum(abs(resUall.^2), 2)/(P-1)/P));          % var mean value over P periods
		varYallm = squeeze((sum(abs(resYall.^2), 2)/(P-1)/P));          % var mean value over P periods
		varYUallm = squeeze(sum(resYall.*conj(resUall), 2)/(P-1)/P);    % covariance mean values over P periods
    end % if P > 1
	Uallm = squeeze(Uallm);
	Yallm = squeeze(Yallm);
	if M == 1
        Uallm = Uallm.';
        Yallm = Yallm.';
        U.mean = Uallm;
        Y.mean = Yallm;
        if P > 1 
            varUallm = varUallm.';
            varYallm = varYallm.';
            varYUallm = varYUallm.';
        end % if P > 1
	end % if M == 1

    if P > 1
        
        % input/output noise (co-)variances
        varUn = mean(varUallm, 1)/M;
        varYn = mean(varYallm, 1)/M;
        varYUn = mean(varYUallm, 1)/M;
        
        % noise covariance matrix [Y; U]
        CYU.n = zeros(2, 2, F);
        CYU.n(1, 1, :) = varYn;
        CYU.n(2, 2, :) = varUn;
        CYU.n(1, 2, :) = varYUn;
        CYU.n(2, 1, :) = conj(varYUn);
        Y.stdn = varYn.^0.5;
        U.stdn = varUn.^0.5;
           
    end % if P > 1

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean and variance over M realisations in the absence of a reference signal          %
    % Averaging spectra over realisations is impossible because they are not synchronized %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin == 2
        
        Gallm = Yallm./Uallm;
        G.mean = mean(Gallm, 1);
        
        % removal noise bias in (cross-)power spectra
        absY2 = mean(abs(Yallm.^2), 1) - M*varYn;
        absU2 = mean(abs(Uallm.^2), 1) - M*varUn;
        YconjU = mean(Yallm .* conj(Uallm), 1) - M*varYUn;
        Y.abs = absY2.^0.5;
        U.abs = absU2.^0.5;
        
        % noise std mean FRF
        G.stdn = abs(G.mean) .* (varYn./absY2 + varUn./absU2 - 2*real(varYUn./YconjU)).^0.5;

		if M > 1
			G.stdNL = std(Gallm, 0, 1) / sqrt(M);       % standard deviation mean value
            Y.stdNL = G.stdNL .* U.abs;                 % standard deviation mean value
		end;
        
    end % if nargin = 2
                         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean and variance over M realisations in the presence of a reference signal %
    % Due to the projection on the reference signal averaging of the spectra over %
    % different realisations is possible                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin == 3
        
        % mean value over realisations
        Y.mean = mean(Yallm, 1);
        U.mean = mean(Uallm, 1);
        Y.abs = abs(Y.mean);
        U.abs = abs(U.mean);
        G.mean = Y.mean./U.mean;
                
        % noise variance mean FRF
        G.stdn = abs(G.mean) .* (varYn./abs(Y.mean.^2) + varUn./abs(U.mean.^2) ...
                                            - 2*real(varYUn./(Y.mean.*conj(U.mean)))).^0.5;
            
       if M > 1
			resU = Uallm - repmat(U.mean, [M, 1]);
			resY = Yallm - repmat(Y.mean, [M, 1]);
			varU = (sum(abs(resU.^2), 1)/(M-1)/M);                  % var mean value over M realisations
			varY = (sum(abs(resY.^2), 1)/(M-1)/M);                  % var mean value over M realisations
			varYU = sum(resY.*conj(resU), 1)/(M-1)/M;               % covariance mean values over M realisations
            
            % total covariance matrix [Y; U]
            CYU.NL = zeros(2, 2, F);
            CYU.NL(1, 1, :) = varY;
            CYU.NL(2, 2, :) = varU;
            CYU.NL(1, 2, :) = varYU;
            CYU.NL(2, 1, :) = conj(varYU);
            Y.stdNL = varY.^0.5;
            U.stdNL = varU.^0.5;
            
            % total std mean FRF
            G.stdNL = abs(G.mean) .* (varY./abs(Y.mean.^2) + varU./abs(U.mean.^2) ...
                                                - 2*real(varYU./(Y.mean.*conj(U.mean)))).^0.5;
        end % if M > 1
        
        
    end % if nargin = 3
    
end % if nargin = 2 or nargin = 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard deviation stochastic nonlinear distortions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (M > 1) & (P > 1)
    
    % estimate GS
    G.stds = (M * abs(G.stdNL.^2 - G.stdn.^2)).^0.5;
    
    % estimate YS
    if nargin > 1   % input/output measurements
        Y.stds = (M * abs(Y.stdNL.^2 - Y.stdn.^2)).^0.5;
    end % if nargin > 1
    
    % estimate US
    if nargin == 3  % reference signal is available
        U.stds = (M * abs(U.stdNL.^2 - U.stdn.^2)).^0.5;
    end % if nargin > 1
    
end



