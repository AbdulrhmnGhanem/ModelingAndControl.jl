function ExerciseDuffing


% Chapter 5 Exercise 79.a and 79.b
% Exercise 79.a: Normal operation, subharmonics, and chaos
%          79.b: Influence initial conditions
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 6 December 2010


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %% 
% Duffing's oscillator                       %%
%                                            %%
%  y'' + y' - 10*y + 100*y^3 = A*cos(3.5*t)  %%
%                                            %%
% Illustration of                            %%
%   1. normal operation   (e.g. A = 0.70)    %%
%   2. period doubling    (e.g. A = 0.81)    %%
%   3. period quadrupling (e.g. A = 0.82)    %%
%   4. chaotic behaviour  (e.g. A = 0.9)     %% 
%                                            %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all;

%
%% The reader has to choose 
%%   -  the excitation        (see lines 46-49)
%%   -  the initial condition (see line 54)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice respone system: select one of the following cases  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Choice = 'period1';         % steady state output has the same period as the input
Choice = 'period2';         % period doubling
Choice = 'period4';         % period quadrupling
Choice = 'chaos';           % chaotic behaviour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the initial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [0.2; -0.2];                       % initial condition position and velocity


%%%%%%%%%%%%%%%%%%%%%
% System parameters %
%%%%%%%%%%%%%%%%%%%%%

theta.a = 1;                % coefficient y'
theta.b = -10;              % coefficient y
theta.c = 100;              % coefficient y^3


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input signal parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Am*cos(Omega*t+Ph)
theta.Omega = 3.5;          % angular frequency input in rad/s; period = 2*pi/Omega
theta.Har = [1];            % fundamental frequency only
theta.Ph = 0;               % zero phase

switch Choice
    case 'period1'
        theta.Am = 0.67;            % the steady state response y(t) has the period 2*pi/Omega
    case 'period2'
        theta.Am = 0.81;            % period doubling: the steady state response y(t)has period 2*(2*pi/Omega)
    case 'period4'
        theta.Am = 0.82;            % period quadrupling: the steady state response y(t)has period 4*(2*pi/Omega)
    case 'chaos'
        theta.Am = 0.90;            % chaotic behaviour
end % switch Choice


%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = theta.Omega/2/pi;      % frequency input           
N = 1024;                   % number of points in one input period
P = 100;                    % number of input periods
fs = N*f0;                  % sampling frequency is chosen such that f0 corresponds to DFT line no. 1
theta.Ts = 1/fs;            % sampling period = initial step size numerical integration
t = [0:1:N*P-1].'/fs;       % time is chosen sufficiently large such that the steady state response
                            % can be observed

                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input signal Duffing's oscillator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = theta.Am*cos(theta.Omega*t + theta.Ph);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution Duffing's oscillator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = SolveDuffing(t, x0, theta);         % numerical integration Duffing's oscillator


%%%%%%%%%%%%%%%%%%%%
% plot the results %
%%%%%%%%%%%%%%%%%%%%

% input and response over the full time interval
figure(1)
subplot(211)
plot(t, u, 'k');
title('Simulated Input and Response', 'FontName','Helvetica LT Std', 'FontSize', 12)
xlabel('Time (s)')
ylabel('Input (a.u.)')
subplot(212);
plot(t, x(:,1), 'k');
TheAxis = axis;
TheAxis(3:4) = [-0.5, 0.5];
axis(TheAxis);
xlabel('Time (s)')
ylabel('Response (a.u.)')
zoom on;
shg

% input and response over the last eight input periods
Plast = 8;
yper = x(end-Plast*N+1:end,1);      % select the last Plast periods
uper = u(end-Plast*N+1:end);        % select the last Plast periods
tper = t(end-Plast*N+1:end);        % select the last Plast periods

figure(2)
subplot(211)
plot(tper, uper, 'k');
title('Last Eight Input Periods of the Simulated Input and Response', ...
     'FontName','Helvetica LT Std', 'FontSize', 11)
xlabel('Time (s)')
ylabel('Input (a.u.)')
subplot(212);
plot(tper, yper, 'k');
TheAxis = axis;
TheAxis(3:4) = [-0.5, 0.5];
axis(TheAxis);
xlabel('Time (s)')
ylabel('Response (a.u.)')
zoom on;
shg

% input and response DFT spectra of the last eight input periods
Y = fft(yper)/N;
Y = Y(1:Plast*N/2);                 % remove the complex conjugate part in the DFT spectrum
U = fft(uper)/N;
U = U(1:Plast*N/2);                 % remove the complex conjugate part in the DFT spectrum
DFTnumber = [1:Plast*N/2].' - 1;    % the selected DFT line numbers
freq = DFTnumber*fs/(Plast*N);      % the corresponding frequencies

figure(3)                           % full DFT spectra
subplot(211)
plot(DFTnumber, db(U),'+');
title('Input and Response DFT Spectra of the Last Eight Input Periods', ...
      'FontName','Helvetica LT Std', 'FontSize', 11)
xlabel('DFT line number')
ylabel('Input (dB)')
subplot(212);
plot(DFTnumber, db(Y),'+');
xlabel('DFT line number')
ylabel('Response (dB)')
zoom on;
shg

figure(4)                           % show first 20 DFT lines
subplot(211)
plot(DFTnumber, db(U),'+');
axis([-0.5 20.5 -350 50]);
title('Zoom of the Input and Response DFT Spectra of the Last Eight Input Periods', ...
      'FontName','Helvetica LT Std', 'FontSize', 11)
xlabel('DFT line number')
ylabel('Input (dB)')
subplot(212);
plot(DFTnumber, db(Y),'+');
axis([-0.5 20.5 -350 50]);
xlabel('DFT line number')
ylabel('Response (dB)')
zoom on;
shg

% phase plane
figure(5)
plot(x(40*N:end,1), x(40*N:end,2), 'k')      % eliminate the transient response
axis([-0.5 0.5 -1 1]);
title('Phase Plane', 'FontName','Helvetica LT Std', 'FontSize', 12)
xlabel('y(t)');
ylabel('dy(t)/dt');
zoom on
shg



function x = SolveDuffing(t, x0, theta);
%
% function x = SolveDuffing(t, x0, theta);
%
%   Solves Duffing's differential equation
%
%       y''(t) + a*y'(t) + b*y(t) + c*y^3(t) = sumk(Ak*cos(k*Omega*t+Fk))
%
%   or in state space form
%
%       x1'(t) = x2(t)
%       x2'(t) = sumk(Ak*cos(k*Omega*t+Fk)) - a*x2(t) - b*x1(t) - c*x1^3(t)
%
%   Output
%       x       =   solution differential equation; size NN x 2
%                       x(:, 1) = y(t)
%                       x(:, 2) = y'(t)
%
%   Input
%       x0      =   initial conditions x0(1) = y(0); x0(2) = y'(0); size 2 x 1
%       theta   =   struct{'a', 'b', 'c', 'Am', 'Ph', 'Omega'}
%                       theta.a     =   a-coefficient Duffing's oscillator
%                       theta.b     =   b-coefficient Duffing's oscillator
%                       theta.c     =   c-coefficient Duffing's oscillator
%                       theta.Am    =   amplitude cosines input; size 1 x numberOfFreq  
%                       theta.Ph    =   phase cosines input; size 1 x numberOfFreq
%                       theta.Har   =   harmonics input; size 1 x numberOfFreq
%                       theta.Omega =   angular frequency cosine input 
%                       theta.Ts    =   initial step size numerical integration 
%      
%
% Rik Pintelon and Johan Schoukens
% January 17, 2007
%

global a b c Am Ph Har Omega Ts

% parameters Duffing's oscillator
a = theta.a;
b = theta.b;
c = theta.c;
Am = theta.Am;
Ph = theta.Ph;
Har = theta.Har;
Omega = theta.Omega;

% set the integration parameters
odeset.RelTol = 1e-10;
odeset.AbsTol = eps;
odeset.InitialStep = Ts/10;
odeset.MaxStep = Ts;

% integration differential equation
[t, x] = ode45(@Duffing, t, x0, odeset);


% ------------------------------------------------------------------------
function xdot = Duffing(t, x);
%
%   function xdot = Duffing(t, x);
%
%           Duffing's oscillator under state space equations form
%               x1'(t) = x2(t)
%               x2'(t) = F*cos(Omega*t) - a*x2(t) - b*x1(t) - c*x1^3(t)
%
%   Output
%       xdot    =	derivative state vector; size 2 x 1
%
%   Input
%       t       =	time
%       x       =	state vector Duffing equation; size 2 x 1
%
%
% Rik Pintelon and Johan Schoukens
% January 16, 2007
%

global a b c Am Ph Har Omega

xdot = zeros(2,1);
xdot(1) = x(2);
xdot(2) = sum(Am.*cos(Har*Omega*t+Ph)) - a*x(2) - b*x(1) - c*x(1)^3;

