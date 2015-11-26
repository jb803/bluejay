%This calculates the equilibrium composition of a mixture undergoing
%multiple dissociation reactions held at constant pressure and temperature

%We must input:
% Initial molar quantities
% Pressure and reference pressure
% A coefficient matrix
% A Kp vector
% A step size vector
%Maximum number of iterations
% Precision threshold (%)


function [equilibN, converged, errorLev] = reactionEquilibrium(N, p, pRef, coeffMat, Kp, epsilon, precThres,maxIters)

    p0= pRef;


    nu = coeffMat;
    
    %We set the number of dissoc. reactions ,k

    k = size(nu,1);

    %We set the number of species, nSpec
    nSpec = length(N);

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

    gamma = Gamma0;

    prec = 100*abs(gamma-Kp)./Kp;


    while max(prec)>precThres && numIters<maxIters
        %We calculate the base coefficients
        nTot = sum(N);
        nTotTerms = nTot*ones(k,1).^(-sumNu);
        NexpTerms = transpose(prod(repmat(N,1,k).^transpose(nu)));
    
        %We evaluate the current Gamma0
        gamma = NexpTerms.*pTerms.*nTotTerms;
        
    
    
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
        
        
        
        %We calculate the base coefficients
        nTot = sum(N);
        nTotTerms = nTot*ones(k,1).^(-sumNu);
        NexpTerms = transpose(prod(repmat(N,1,k).^transpose(nu)));
    
        %We evaluate the current Gamma0
        gamma = NexpTerms.*pTerms.*nTotTerms
        
    
        prec = 100*abs(gamma-Kp)./Kp;
    
        numIters = numIters+1;
        
        errorLev = prec;
            
    
    
    end
    numIters
    
    if max(prec)<precThres 
        converged=true;
    else
        converged = false;
    end
    
    equilibN = N;

end


