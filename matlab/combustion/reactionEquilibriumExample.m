%This calculates the equilibrium composition of a mixture undergoing
%multiple dissociation reactions held at constant pressure and temperature

%For this the species are:
%1 -> CO
%2 -> O2
%3 -> CO2
%4 -> H20
%5 -> H2
%6 -> N2
%7 -> N
%8 -> OH
%9 -> NH3
%10-> NO

%This undergoes 2 dissoc. reactions
%-CO -0.5O2 + CO2 = 0 R7
%-CO - H20  + CO2 + H2 = 0 R8

tic

%Define pressure and ref pressure
%p = 12; %bar
%p0= 1; %bar

%We set up the initial molar quantites
%N = [ 1; 3; 5; 5; 1; 2; 0.1; 0.01; 0.01; 1];

%We set up the coefficient matrix. Each row follows a seperate reactant
%nu = [ -1  -0.5 1  0 0    0    0  0  0 0 %R7
%       -1    0  1 -1 1    0    0  0  0 0%R8
%        0    0  0  0 0    1    -2 0  0 0%R2
%        0    0  0  1 -0.5 0    0 -1 0  0%R6
%        0    0  0  0 -1.5 -0.5 0  0 1  0%R9
%        0    1  0  0 0    1    0  0 0  -2];%R4
%We set up the Kp vector for T=2000K

%Kp = [exp(6.634); exp(-1.510); exp(41.645);exp(8.727); exp(-10.810); exp(7.824)];

%We set up the step sizes
%epsilon = [ 1/1000000; 1/10000; 1E-21;1/1000000; 1/10; 1/1000000 ];




%2nd example
%1->N2
%2->CO2
%3->CO
%4->O2

N       = [1; 1;   1]; %Initial species vector
nu      = [1 -1  -0.5]  %reqction 7
Kp = exp(6.634);
p0 = 1;
p = 2;
epsilon = 1/1000000;

[N converged errorLev] = reactionEquilibrium(N,p,p0,nu,Kp,epsilon,0.05,30000)


toc


