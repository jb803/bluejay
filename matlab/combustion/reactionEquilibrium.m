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
p = 12; %bar
p0= 1; %bar

%We set up the initial molar quantites
N = [ 1; 3; 5; 5; 1; 2; 0.1; 0.01; 0.01; 1];
%N = [ 2; 2; 2; 2; 2; 2; 2; 2; 2; 2];
%We set up the coefficient matrix. Each row follows a seperate reactant
nu = [ -1  -0.5 1  0 0    0    0  0  0 0 %R7
       -1    0  1 -1 1    0    0  0  0 0%R8
        0    0  0  0 0    1    -2 0  0 0%R2
        0    0  0  1 -0.5 0    0 -1 0  0%R6
        0    0  0  0 -1.5 -0.5 0  0 1  0%R9
        0    1  0  0 0    1    0  0 0  -2];%R4
%We set up the Kp vector for T=2000K

Kp = [exp(6.634); exp(-1.510); exp(41.645);exp(8.727); exp(-10.810); exp(7.824)];

%We set up the step sizes
epsilon = [ 1/1000000; 1/10000; 1E-21;1/1000000; 1/10; 1/1000000 ];

%We set the number of dissoc. reactions ,k

k = 6;

%We set the number of species, nSpec
nSpec = 10;

%We prepare the compliance matrix C
C=zeros(nSpec,k);


%We set up the balance matrix
for j=1:k
    B = zeros(nSpec-1,nSpec);
    iStarFound = false;
    counter = 1;
    iStar=0;
    for i=1:nSpec
        if nu(j,i) == 0
            B(counter,i)=1;
            counter=counter+1;
        elseif iStarFound == false
            iStarFound=true;
            iStar = i;
        else
            B(counter,i)=-1/nu(j,i);
            B(counter,iStar)=1/nu(j,iStar);
            counter=counter+1;
        end
        
    end
    C(:,j)=null(B,'r');
end

numIters=0;

%We calculate the coefficient sums
sumNu=transpose(sum(transpose(nu)));

pTerms = ((p/p0)*ones(k,1)).^sumNu;

%We calculate the initial gammas
nTot = sum(N);
nTotTerms = nTot*ones(k,1).^(-sumNu);
NexpTerms = transpose(prod(repmat(N,1,k).^transpose(nu)));
Gamma0=NexpTerms.*pTerms.*nTotTerms;

maxStepSize = epsilon.*(Kp-Gamma0);

gHist = Gamma0;

precThres=0.05;

gamma = Gamma0;

prec = 100*abs(gamma-Kp)./Kp;


while max(prec)>precThres
    %We calculate the base coefficients
    nTot = sum(N);
    nTotTerms = nTot*ones(k,1).^(-sumNu);
    NexpTerms = transpose(prod(repmat(N,1,k).^transpose(nu)));
    
    %We evaluate the current Gamma0
    gamma = NexpTerms.*pTerms.*nTotTerms;
    
    gHist=[gHist gamma];
    
    %We now need to generate the equilbrium matrix M
    xi=sumNu/nTot;
    
    E = zeros(k,nSpec);
    
    for j=1:k
        for i=1:nSpec
            if nu(j,i)== 0
                E(j,i) = 0;
            else
                E(j,i)=(nu(j,i)/N(i))-xi(j);
            end
        end
    end
    G = zeros(k,1);
    %We calculate G the step matrix
    for j = 1:k
        if abs(Kp(j)-gamma(j))>maxStepSize(j)
            G(j) = abs(maxStepSize(j)) * sign(Kp(j)-gamma(j));
        else
            G(j) = Kp(j) - gamma(j);
        end
    end
            
    
    reactionDelta = (E*C)\G;
    
    speciesDelta = C*reactionDelta;
    
    N = N + speciesDelta;
    
    prec = 100*abs(gamma-Kp)./Kp;
    
    numIters = numIters+1;
            
    
    
end

%plot(gHist(1,:));
%title('R7');
%figure
%plot(gHist(2,:));
%title('R8');
%figure
%plot(gHist(3,:));
%title('R2');
%figure
%plot(gHist(4,:));
%title('R6');
%figure
%plot(gHist(5,:));
%title('R9');
%figure
%plot(gHist(6,:));
%title('R4');

toc

numIters


