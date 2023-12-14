% Chapter 5 Exercise 78
% Uniform versus pointwise convergence
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 6 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform versus non-uniform polynomial approximation %%
%                                                     %%
% Illustrated on                                      %%
%   1. atan function                                  %%
%   2. step function                                  %%
%                                                     %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% atan function                               %
%                                             %
%   1. Taylor series expansion                %
%   2. Polynomial Least squares approximation %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 400;                                                % number of data points
DegreeAll = [1, 3, 7, 15, 23];                          % odd degrees of the polynomial approximation
nDeg = length(DegreeAll);                               % number of degrees considered
u = linspace(-3, 3, N).';                               % the points are linearly spaced in [-3, 3]
yatan = atan(u);                                        % the exact value

% Taylor series approximation
ytaylor = zeros(N, nDeg);
CoeffTaylor = zeros(1, max(DegreeAll)+1);
kk = [1:2:DegreeAll(end)];                              % degrees of all the monomials in the highest degree expansion
CoeffTaylor(2:2:end) = ((-1).^((kk-1)/2)) ./ kk;        % coefficients of odd degree monomials (even degree are zero)
for ii = 1:nDeg
    Coeff = fliplr(CoeffTaylor(1:DegreeAll(ii)+1));     % coefficients taylor approximation of degree DegreeAll(ii)
    ytaylor(:, ii) = polyval(Coeff, u);                 % value of the Taylor series expansion
end % ii
eTaylor = ytaylor - repmat(yatan, [1, nDeg]);           % approximation error of the Taylor series expansion

% Polynomial least squares approximation
yls = zeros(N, nDeg);
uls = u/std(u, 1);                                      % normalisation to improve the numerical conditioning
for ii = 1:nDeg
    Coeffls = polyfit(uls, yatan, DegreeAll(ii));       % coefficients least squares approximation of degree DegreeAll(ii)
    yls(:, ii) = polyval(Coeffls, uls);                 % value of the least squares approximation
end % ii
els = yls - repmat(yatan, [1, nDeg]);                   % approximation error of the least squares fit

% Plot the results
figure(1)
subplot(221)
plot(u, ytaylor, 'r', u, yatan,'k')
axis([-4, 4, -2, 2])
xlabel('u');
ylabel('y');
title('Taylor Approximation atan', 'FontName','Helvetica LT Std', 'FontSize', 10);

subplot(222)
plot(u, yls, 'r', u, yatan,'k')
xlabel('u');
ylabel('y');
title('Least Squares Approximation atan', 'FontName','Helvetica LT Std', 'FontSize', 10);

subplot(223)
plot(u, db(eTaylor), 'r')
% axis([-4, 4, 1e-20, 1e20])
xlabel('u');
ylabel('error (dB)');

subplot(224)
plot(u, db(els), 'r')
% axis([-4, 4, 1e-20, 1e20])
xlabel('u');
ylabel('error (dB)');


%%%%%%%%%%%%%%%%%%%%%%%%
% triangular function  %
%%%%%%%%%%%%%%%%%%%%%%%%

N = 400*3;                                              % number of data points
DegreeAll = [2, 6, 40];                                 % odd degrees of the polynomial approximation
nDeg = length(DegreeAll);                               % number of degrees considered
ut = linspace(-3, 3, N).';                              % the points are linearly spaced in [-3, 3]
utn = ut/std(ut, 1);                                    % normalised input (to improve the numerical conditioning)
yt = abs(ut) - 1.5;                                     % the exact value

% Polynomial least squares approximation
ytls = zeros(N, nDeg);
for ii = 1:nDeg
    Coeffls = polyfit(utn, yt, DegreeAll(ii));          % coefficients least squares approximation of degree DegreeAll(ii)
    ytls(:, ii) = polyval(Coeffls, utn);                % value of the least squares approximation
end % ii
etls = ytls - repmat(yt, [1, nDeg]);                    % approximation error of the least squares fit

% Plot the results of the least squares approximation
figure(2)
subplot(211)
plot(ut, ytls, 'r', ut, yt,'k')
axis([-3, 3,-2, 2]); 
xlabel('u');
ylabel('y');
title('Least Squares Approximation Triangle', 'FontName','Helvetica LT Std', 'FontSize', 12);

subplot(212)
plot(ut, db(etls), 'r')
axis([-3, 3,-100, 0]); 
xlabel('u');
ylabel('error (dB)');


%%%%%%%%%%%%%%%%%%
% step function  %
%%%%%%%%%%%%%%%%%%

N = 400*3;                                              % number of data points
DegreeAll = [1, 7, 41];                                 % odd degrees of the polynomial approximation
nDeg = length(DegreeAll);                               % number of degrees considered
ud = linspace(-3, 3, N).';                              % the points are linearly spaced in [-3, 3]
udn = ud/std(ud, 1);                                    % normalised input (to improve the numerical conditioning)
yd = sign(ud);                                          % the exact value

% Polynomial least squares approximation
ydls = zeros(N, nDeg);
for ii = 1:nDeg
    Coeffls = polyfit(udn, yd, DegreeAll(ii));          % coefficients least squares approximation of degree DegreeAll(ii)
    ydls(:, ii) = polyval(Coeffls, udn);                % value of the least squares approximation
end % ii
edls = ydls - repmat(yd, [1, nDeg]);                    % approximation error of the least squares fit

% Plot the results of the least squares approximation
figure(3)
subplot(211)
plot(ud, ydls, 'r', ud, yd,'k')
axis([-3, 3,-2, 2]); 
xlabel('u');
ylabel('y');
title('Least Squares Approximation Step', 'FontName','Helvetica LT Std', 'FontSize', 12);

subplot(212)
plot(ud, db(edls), 'r')
xlabel('u');
ylabel('error (dB)');



