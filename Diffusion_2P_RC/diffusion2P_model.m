%% A simple diffusion+ECN model

function [voltOut]=diffusion2P_model(k,currData,timeData,socData,tempData,OcvLuts)
%params

N=1; % controls timestep
dt=1/N;
Tref=273+25;

SoC0=socData(1);
voltOut=ones(length(timeData),1);
Vrc=0;

R_0=k(1);
R_1 = k(2);
C_1=k(3);
tauP_0=k(4);
tauN_0=k(5);
kd=k(6);
Ea1=k(7);
Ea2=k(8);
Ea3=k(9);


% diffusion settings
Q=4.7455*3600; %capacity, should be an argument but lazy for the moment
R = 1; % particle radius [m]
Nr = 8; % number of "shells" radially
dR = R/Nr; % width of each "shell"
Sa = 4*pi*(R*(1:Nr)/Nr).^2; % outer surface area of each shell
dV = (4/3)*pi*((R*(1:Nr)/Nr).^3-(R*(0:Nr-1)/Nr).^3); % vol. of ea. shell
SoCN = SoC0*ones(1,(Nr)); % 'concentration' profile versus "r" dimension for neg electrode
SoCP=(1-SoC0)*ones(1,Nr);
SoCNs = zeros(size(timeData)); % concentration at surface
SoCPs=zeros(size(timeData)); 
% SoCavg=SoC0*ones(size(timeData));
SoCNs(1) = SoC0;
SoCPs(1)=1-SoC0;

h(1)=-1;
k_hyst=10;
hyst=OcvLuts.Components.hystAmp(:,5); %Hyst data parameters, fifth column for 25 deg.
hyst_0=OcvLuts.Components.hystInst(:,5);
% SoCr=ones(length(timeData),Nr); %internal SoC, maybe useful for gradient]

%if using finer timestep 'interploate' input current array to larger size
times=N*length(timeData);
curr_Data=currData;

% for ii=2:length(currData)
%   curr_Data(N*ii+1:N*ii+10)=currData(ii);
% end
% curr_Data(1:N)=currData(1);

% calc diffusion 
for timestep = 1:times
  
    R0=R_0*exp(-Ea1/8.314*(-1/(273+tempData(timestep))+1/Tref));
    R1=R_1*exp(-Ea2/8.314*(-1/(273+tempData(timestep))+1/Tref));
    tauP=tauP_0*exp(-Ea3/8.314*(-1/(273+tempData(timestep))+1/Tref));
    tauN=tauN_0*exp(-Ea3/8.314*(-1/(273+tempData(timestep))+1/Tref));
    
    IR0=R0.*curr_Data(timestep);     
    Vrc=R1*(exp(-dt/R1/C_1).*(Vrc/R1)+(1-exp(-dt/R1/C_1)).*curr_Data(timestep));

    M_hyst=interp1(OcvLuts.Dims.soc,hyst,socData(timestep));
    M0=interp1(OcvLuts.Dims.soc,hyst_0,socData(timestep));
    h=exp(-dt*abs(curr_Data(timestep))*k_hyst/(Q))*h+sign(curr_Data(timestep))*(1-exp(-dt*abs(curr_Data(timestep)*k_hyst/(Q))));
    U_hyst=M_hyst.*h+sign(curr_Data(timestep)).*M0;

    % tau=tau0/(socData(timestep)+beta); 
    % tau=beta*(1-socData(timestep))+tau0; %quadratic form
    fluxN = -1/tauN*diff(SoCN)/dR; % flux at surfaces between "bins"
    fluxP=  -1/tauP*diff(SoC)/dR;
    MN = flux.*Sa(1:end-1); % total SoC crossing surface between bins
    MP = flux.*Sa(1:end-1); 
    SoCN= SoCN+ ([0 MN] - [MN 0])*dt./dV; % conc. change via diffusion
    SoCP= SoCP+ ([0 MP] - [MP 0])*dt./dV; % conc. change via diffusion
    % SoCr(timestep,:)=SoC;
    SoCN(end) = SoCN(end) + (curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end); % at boundary neg elec
    SoCP(end) = SoCP(end) -(curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end);  % at boundary pos elec
    SoCNs(timestep) = min(1,SoCN(end)); % surface soc
    SoCPs(timestep) = min(1,SoCP(end));
    % SoCavg(timestep)= socData(timestep);
    % OCVcell(timestep)=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.ocv,SoCavg(timestep),tempData(timestep),'makima');% ocv at avearge soc
    OCVP_surf(timestep)=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.U_p,1-SoCs(timestep),tempData(timestep),'makima'); %ocv at surface soc
    OCVN_surf(timestep)=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.U_n,SoCs(timestep),tempData(timestep),'makima');
    %Vdiff(timestep)=-(kd*(OCVcell(timestep)-OCVcell_surf(timestep)));
    %v_Out(timestep)=IR0+Vrc+Vdiff(timestep)+OCVcell(timestep)+U_hyst;
    v_Out(timestep)=IR0+Vrc+U_hyst+kd*(OCVP_surf-OCVN_surf);
    end

    voltOut=v_Out';

    % for i=2:length(timeData)
    %     
    % voltOut(i)=v_Out(N*(i-1)+1);
    % end
    % voltOut(1)=v_Out(1);


end

