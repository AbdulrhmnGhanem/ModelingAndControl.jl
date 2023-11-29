% Chapter 5 Exercise 81
% Detection, quantification, and classification of nonlinearities
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %% 
% Detection, qualification, and classification nonlinearities  %%
% using a odd multisine where some of the odd harmonics are    %%
% not excited                                                  %%
%                                                              %%
%   System considered                                          %%
%                                                              %%
%  y''/w0^2 + 2*zeta/w0*y' + y = u + alpha*u'^2 + beta*u^2*u'  %%
%                                                              %%
%   Input                                                      %%
%          1. odd multisine where one of the odd               %%
%             harmonics is not excited                         %%
%          2. odd-odd multisine where every second odd         %%
%             harmonic is not excited                          %%
%                                                              %%
% Rik Pintelon and Johan Schoukens                             %% 
% March 2, 2007                                                %%
%                                                              %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition left hand side NL dynamic system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = 30;                                % resonance angular frequency
zeta = 0.2;                             % damping ratio
alpha = 1e-7;
beta = 1e-3;


%%%%%%%%%%%%%%%%%%%%%%%%%
% multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%

% basic parameters multisine
fo = 1;                                 % frequency odd multisine
foo = 1/2;                              % frequency odd-odd multisine (double frequency resolution)
fs = 200;                               % sampling frequency
No = fs/fo;                             % number of points in one period of the odd multisine
Noo = fs/foo;                           % number of points in one period of the odd-odd multisine
to = [0:No-1].';                        % time in samples of the odd multisine
too = [0:Noo-1].';                      % time in samples of the odd-odd multisine

% odd multisine with 4 components: [1, 3, 5, 7] * f0
Uo = zeros(No, 1) + eps;                % add eps to make the non-excited lines visible in a dB-scale
kko = [1; 3; 5; 7];
Uo(kko+1) = exp(-sqrt(-1)*pi*kko.*(kko-1)/length(kko));
uo = real(ifft(Uo));
uo = uo/std(uo, 1);                     % rms value equal to one

% odd-odd multisine with 4 components: [1, 5, 9, 13] * f0
Uoo = zeros(Noo, 1) + eps;              % add eps to make the non-excited lines visible in a dB-scale
kkoo = [1; 5; 9; 13];
Uoo(kkoo+1) = exp(-sqrt(-1)*pi*kkoo.*(kkoo-1)/length(kkoo));
uoo = real(ifft(Uoo));
uoo = uoo/std(uoo, 1);                  % rms value equal to one


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivative odd multisine
Uo = fft(uo);
dUo = zeros(No, 1);
dUo(kko+1) = Uo(kko+1).*(2*pi*fo*sqrt(-1)*kko);
duo = 2*real(ifft(dUo));

% derivative odd-odd multisine
Uoo = fft(uoo);
dUoo = zeros(Noo, 1);
dUoo(kkoo+1) = Uoo(kkoo+1).*(2*pi*foo*sqrt(-1)*kkoo);
duoo = 2*real(ifft(dUoo));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response of the nonlinear dynamic system to odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% right hand side differential equation for odd multisine: z = u + alpha*(u'^2 + u^2*u')
zo = uo + alpha*duo.^2 + beta*uo.^2 .* duo;       % odd multisine
Zo = fft(zo);

% DFT lines at which output is calculated (frequencies do not exceed 3*f0)
ko = [0:1:25].';

% transfer function
Go = 1 ./ polyval([1/w0^2, 2*zeta/w0, 1], 2*pi*fo*sqrt(-1)*ko);

% output spectrum
Yo = zeros(No, 1);
Yo(ko+1) = Go .* Zo(ko+1);
yo = 2*real(ifft(Yo));

% selection DFT lines shown in the figures
DFTlineso = [0:28].';
freqo = DFTlineso*fs/No;                     % frequencies odd-odd multisine
Uo = Uo(DFTlineso+1)/No + eps;               % add eps to show the non-excited lines in a dB scale
Yo = Yo(DFTlineso+1)/No + eps;

% selection even, odd excited, and odd non-excited harmonics
Harmo = struct('odd',[], 'even', []);        % odd, even	=> respectively odd and even harmonics
Harmo.odd = struct('E', [], 'NE', []);       % E, NE     => respectively excited and non-excited harmonics
Harmo.even = DFTlineso(1:2:end);
Harmo.odd.E = kko;
Harmo.odd.NE = DFTlineso(2:2:end);
Fexco = length(Harmo.odd.E);
RemoveIndexo = zeros(1, Fexco);
for ii = 1:Fexco
    RemoveIndexo(ii) = find(Harmo.odd.NE == Harmo.odd.E(ii));
end % ii
Harmo.odd.NE(RemoveIndexo) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response of the nonlinear dynamic system to odd-odd multisine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% right hand side differential equation for odd multisine: z = u + alpha*(u'^2 + u^2*u')
zoo = uoo + alpha*duoo.^2 + beta*uoo.^2 .* duoo;       % odd-odd multisine
Zoo = fft(zoo);

% DFT lines at which output is calculated (frequencies do not exceed 3*f0)
koo = [0:1:50].';

% transfer function
Goo = 1 ./ polyval([1/w0^2, 2*zeta/w0, 1], 2*pi*foo*sqrt(-1)*koo);

% output spectrum
Yoo = zeros(Noo, 1);
Yoo(koo+1) = Goo .* Zoo(koo+1);
yoo = 2*real(ifft(Yoo));

% selection DFT lines shown in the figures
DFTlinesoo = [0:56].';
freqoo = DFTlinesoo*fs/Noo;                    % frequencies odd-odd multisine
Uoo = Uoo(DFTlinesoo+1)/Noo + eps;             % add eps to show the non-excited lines in a dB scale
Yoo = Yoo(DFTlinesoo+1)/Noo + eps;

% selection even, odd excited, and odd non-excited harmonics
Harmoo = struct('odd',[], 'even', []);        % odd, even	=> respectively odd and even harmonics
Harmoo.odd = struct('E', [], 'NE', []);       % E, NE     => respectively excited and non-excited harmonics
Harmoo.even = DFTlinesoo(1:2:end);
Harmoo.odd.E = kkoo;
Harmoo.odd.NE = DFTlinesoo(2:2:end);
Fexcoo = length(Harmoo.odd.E);
RemoveIndexoo = zeros(1, Fexcoo);
for ii = 1:Fexcoo
    RemoveIndexoo(ii) = find(Harmoo.odd.NE == Harmoo.odd.E(ii));
end % ii
Harmoo.odd.NE(RemoveIndexoo) = [];


%%%%%%%%%%%%%%%%%%%%
% Plot the results %
%%%%%%%%%%%%%%%%%%%%

% odd multisine
figure(1)
subplot(221)
plot(to, uo)
axis([-10, 210, -3, 3]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(to, yo)
axis([-10, 210, -5.5, 5.5]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(freqo, db(Uo), 'k+')
axis([-1, 32, -340, 40]);
xlabel('Frequency (Hz)')
ylabel('Input (dB)')

subplot(224)
plot(freqo(Harmo.odd.E+1), db(Yo(Harmo.odd.E+1)), 'k+', freqo(Harmo.odd.NE+1), db(Yo(Harmo.odd.NE+1)), 'ro', ...
     freqo(Harmo.even+1), db(Yo(Harmo.even+1)), 'g*');
axis([-1, 32, -340, 40]);
xlabel('Frequency (Hz)')
ylabel('Output (dB)')

annotation('textbox',[0.34 0.8 0.05 0.2],'Color','k','String','Response to Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

% odd-odd multisine
figure(2)
subplot(221)
plot(too, uoo)
axis([-10, 210, -3, 3]);
xlabel('Samples')
ylabel('Input (a.u.)')

subplot(222)
plot(too, yoo)
axis([-10, 210, -5.5, 5.5]);
xlabel('Samples')
ylabel('Output (a.u.)')

subplot(223)
plot(freqoo, db(Uoo), 'k+')
axis([-1, 32, -340, 40]);
xlabel('Frequency (Hz)')
ylabel('Input (dB)')

subplot(224)
plot(freqoo(Harmoo.odd.E+1), db(Yoo(Harmoo.odd.E+1)), 'k+', freqoo(Harmoo.odd.NE+1), db(Yoo(Harmoo.odd.NE+1)), 'ro', ...
     freqoo(Harmoo.even+1), db(Yoo(Harmoo.even+1)), 'g*');
axis([-1, 32, -340, 40]);
xlabel('Frequency (Hz)')
ylabel('Output (dB)')

annotation('textbox',[0.3 0.8 0.05 0.2],'Color','k','String','Response to Odd-Odd Multisine',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

