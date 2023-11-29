
% Chapter 1 Exercise 5.b
% Weighted least squares estimation: A study of the variance
%
% This excercise starts from the data generated in Exercise 6
%
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 24 November 2010

disp('Generate first the raw data by running Exercise 6')


varLS=((i'*i)^-1) *(i'*(i.*w))*(i'*i)^-1;   % calculate the theoretical variance of the LS
varWLS=(i'*(i./w))^-1;                      % calculate the theoretical variance of the WLS

% display the results
'theoretical variance LS',disp(sqrt(varLS))
'experimental variance LS',disp(std(R1))
'theoretical variance WLS',disp(sqrt(varWLS))
'experimental variance WLS',disp(std(R2))


