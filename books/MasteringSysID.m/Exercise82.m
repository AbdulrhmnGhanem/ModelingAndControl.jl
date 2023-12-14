% Chapter 5 Exercise 82
% Influence DC values signals on the linear approximation
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 7 December 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influence of the DC values of the input/output      %%
% signals on the linear approximation of a static     %%
% nonlinear system.                                   %%
%                                                     %%
% Illustrated on y(t) = u^3(t) where                  %%
%   1. u(t) is uniform                                %%
%   2. u(t) is gaussian                               %%
%                                                     %%
% Rik Pintelon and Johan Schoukens                    %% 
% March 28, 2007                                      %%
%                                                     %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% Input/output signals %
%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian input signal
N = 40000;                                                          % number of data points
mu = 1.5;                                                           % mean value
sigma = 1/sqrt(12);                                                 % standard deviation
uu = rand(N, 1) + 1;                                                % white uniform noise
ug = randn(N, 1)*sigma + mu;                                        % white Gaussian noise
uDC_u = mean(uu);
uDC_g = mean(ug);

% output signals
yu = uu.^3;                                                         % output signal used for the LA's
yg = ug.^3;                                                         % output signal used for the LA's
yDC_u = mean(yu);
yDC_g = mean(yg);

% input/output signals for plotting the true nonlinear function
u = linspace(1, 2, N).';
y = u.^3;                                   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA without DC values removed in input/output signals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coeff_DC_u = sum(yu.*uu)./sum(uu.^2);                               % coefficients best linear approximation (BLA)
BLA_DC_u = Coeff_DC_u*u;                                            % value of the BLA
Coeff_DC_g = sum(yg.*ug)./sum(ug.^2);                               % coefficients best linear approximation (BLA)
BLA_DC_g = Coeff_DC_g*u;                                            % value of the BLA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLA without DC values removed in input/output signals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coeff_noDC_u = sum((yu-yDC_u).*(uu-uDC_u))./sum((uu-uDC_u).^2);     % coefficients best linear approximation (BLA)
BLA_noDC_u = Coeff_noDC_u*(u-uDC_u) + yDC_u;                        % value of the BLA
Coeff_noDC_g = sum((yg-yDC_g).*(ug-uDC_g))./sum((ug-uDC_g).^2);     % coefficients best linear approximation (BLA)
BLA_noDC_g = Coeff_noDC_g*(u-uDC_g) + yDC_g;                        % value of the BLA


%%%%%%%%%%%%%%%%%%%%
% Plot the results %
%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf, 'Position',[50 350 1200 400])
subplot(121)
plot(u, y, 'k', u, BLA_DC_u, 'r', u, BLA_noDC_u, 'b')
% axis([-4, 4, -2, 2])
xlabel('u(t)');
ylabel('y(t)');

subplot(122)
plot(u, y, 'k', u, BLA_DC_g, 'r', u, BLA_noDC_g, 'b')
% axis([-4, 4, -2, 2])
xlabel('u(t)');
ylabel('y(t)');

annotation('textbox',[0.25 0.8 0.05 0.2],'Color','k','String','Linear approximations without (red) and with (blue) DC removal in the input/output signals',...
           'FontName','Helvetica LT Std','FontSize', 12 , 'LineStyle', 'none','FitHeightToText','on',...
           'HorizontalAlignment', 'center' );

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison estimates - asymptotic values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coeff0_noDC_g = 3*sigma^2 + 3*mu^2;
Coeff0_noDC_u = 9/5*sigma^2 + 3*mu^2;
Coeff0_DC_g = (3*sigma^4 + mu^4 + 6*mu^2*sigma^2)/(sigma^2 + mu^2);
Coeff0_DC_u = (9/5*sigma^4 + mu^4 + 6*mu^2*sigma^2)/(sigma^2 + mu^2);

fprintf('\nLinear approximation for uniform input (estimate, true value): %e %e \n', Coeff_DC_u, Coeff0_DC_u);
fprintf('Linear approximation for Gaussian input (estimate, true value): %e %e \n', Coeff_DC_g, Coeff0_DC_g);
fprintf('Best linear approximation for uniform input (estimate, true value): %e %e \n', Coeff_noDC_u, Coeff0_noDC_u);
fprintf('Best linear approximation for Gaussian input (estimate, true value): %e %e \n', Coeff_noDC_g, Coeff0_noDC_g);
