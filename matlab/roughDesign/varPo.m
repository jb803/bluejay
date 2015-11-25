%Gas properties

gamma =1.333;
cp=1300;
T0=1800;

%Set up property table
%M pPo VcpTo mDotCpToOnAPo
propTable=[];
for M=0.01:0.01:6
    pPo=(1+0.5*(gamma-1)*M^2)^(-gamma/(gamma-1));
    VcpTo=sqrt(gamma-1)*M*(1+0.5*(gamma-1)*M^2)^-0.5;
    mDotCpTo=(gamma/sqrt(gamma-1))*M*(1+0.5*(gamma-1)*M^2)^(-0.5*((gamma+1)/(gamma-1)));
    
    propTable = [propTable;
                 M pPo VcpTo mDotCpTo];
end

F1 = (gamma/sqrt(gamma-1))*1*(1+0.5*(gamma-1)*1)^(-0.5*((gamma+1)/(gamma-1)));
    
%We now do p/p0
desiredImpulse = 7450; %NS
Thrust = 1000;
percMass=zeros(length(Thrust),2);
for j=1:length(Thrust)
    pStag = [10 20];
    ARatio =zeros(1,length(pStag));
    ReqMDot = zeros(1,length(pStag));
    MaxVehicleMass=Thrust/50;
    for i=1:length(pStag)
        pRatio=1/pStag(i);
        M=interp1(propTable(:,2),propTable(:,1),pRatio);
        V=sqrt(cp*T0)*interp1(propTable(:,1),propTable(:,3),M);
        F=interp1(propTable(:,1),propTable(:,4),M);
        ARatio(i)=F1/F;
        ReqMDot(i)=Thrust(j)/V;
        ReqBurnTime=desiredImpulse/Thrust(j);
        ReqMass=ReqMDot(i)*ReqBurnTime
        percMass(j,i)=100*ReqMass/MaxVehicleMass(j);
        
 
    end
end
figure
plot(Thrust,percMass(:,1),'b-',Thrust,percMass(:,2),'r-');
xlabel('Thrust /N');
ylabel('Fuel mass as percentage maximum vehicle mass (%)');
legend('10bar','20 bar');
title('Fuel requirements as a percentage of maximum vehicle mass based on TWR')
figure
[ax,line1,line2]=plotyy(pStag,ReqMDot,pStag,ARatio);
xlabel('Chamber Stagnation Pressure (bar)');
ylabel(ax(1),'Required $\dot{m}$','Interpreter','LaTex');
ylabel(ax(2),'Area Ratio');
title('Chamber pressure variation');