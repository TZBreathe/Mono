%% A simple diffusion+ECN model

function [voltOut]=diffusion2P_model_run_lut(k,currData,timeData,socData,tempData,OcvLuts)
%params

N=1; % controls timestep
dt=1;
Tref=273+25;
C_lut=[0.1; 0.5; 1; 2; 3; 4; 5].*4.8;
%C_lut=[12.5; 15; 20];
SoC0=socData(1);
%v_Out=ones(length(timeData),1)';
%voltOut=ones(length(timeData),1);
Vrc=0;

% diffusion settings
Q=4.7455*3600; %capacity, should be an argument but lazy for the moment
R = 1; % particle radius [m]
Nr = 8; % number of "shells" radially
dR = R/Nr; % width of each "shell"
Sa = 4*pi*(R*(1:Nr)/Nr).^2; % outer surface area of each shell
dV = (4/3)*pi*((R*(1:Nr)/Nr).^3-(R*(0:Nr-1)/Nr).^3); % vol. of ea. shell
SoCN = SoC0*ones(1,(Nr)); % 'concentration' profile versus "r" dimension for neg electrode
SoCP=(1-SoC0)*ones(1,Nr);
% SoCavg=SoC0*ones(size(timeData));
SoCNs = SoC0;
SoCPs = 1-SoC0;

h(1)=-1;
k_hyst=10;
hyst=OcvLuts{1}.Components.hystAmp(:,5); %Hyst data parameters, fifth column for 25 deg.
hyst_0=OcvLuts{1}.Components.hystInst(:,5);
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
     
  R_0=interp1(C_lut,k(1,:),curr_Data(timestep),'linear','extrap');
  R_1 = interp1(C_lut,k(2,:),curr_Data(timestep),'linear','extrap');
  C_1=interp1(C_lut,k(3,:),curr_Data(timestep),'linear','extrap');
  tauP_0=interp1(C_lut,k(4,:),curr_Data(timestep),'linear','extrap');
  tauN_0=interp1(C_lut,k(5,:),curr_Data(timestep),'linear','extrap');
  kd=interp1(C_lut,k(6,:),curr_Data(timestep),'linear','extrap');
  Ea1=k(7,1);
  Ea2=k(8,1);
  
    
    
    
    R0=R_0*exp(-Ea1/8.314*(-1/(273+tempData(timestep))+1/Tref));
    R1=R_1*exp(-Ea1/8.314*(-1/(273+tempData(timestep))+1/Tref));
    tauP=tauP_0*exp(-Ea2/8.314*(-1/(273+tempData(timestep))+1/Tref));
    tauN=tauN_0*exp(-Ea2/8.314*(-1/(273+tempData(timestep))+1/Tref));
    
    IR0=R0.*curr_Data(timestep);     
    Vrc=R1*(exp(-dt/R1/C_1).*(Vrc/R1)+(1-exp(-dt/R1/C_1)).*curr_Data(timestep));

    M_hyst=interp1(OcvLuts{1}.Dims.soc,hyst,socData(timestep));
    M0=interp1(OcvLuts{1}.Dims.soc,hyst_0,socData(timestep));
    h=exp(-dt*abs(curr_Data(timestep))*k_hyst/(Q))*h+sign(curr_Data(timestep))*(1-exp(-dt*abs(curr_Data(timestep)*k_hyst/(Q))));
    U_hyst=M_hyst.*h+sign(curr_Data(timestep)).*M0;

    % tau=tau0/(socData(timestep)+beta); 
    % tau=beta*(1-socData(timestep))+tau0; %quadratic form
    fluxN = -1/tauN*diff(SoCN)/dR; % flux at surfaces between "bins"
    fluxP=  -1/tauP*diff(SoCP)/dR;
    MN = fluxN.*Sa(1:end-1); % total SoC crossing surface between bins
    MP = fluxP.*Sa(1:end-1); 
    SoCN= SoCN+ ([0 MN] - [MN 0])*dt./dV; % conc. change via diffusion
    SoCP= SoCP+ ([0 MP] - [MP 0])*dt./dV; % conc. change via diffusion
    % SoCr(timestep,:)=SoC;
    SoCN(end) = SoCN(end) + (curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end); % at boundary neg elec
    SoCP(end) = SoCP(end) -(curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end);  % at boundary pos elec
    SoCNs = min(1,SoCN(end)); % surface soc
    SoCPs = min(1,SoCP(end));
    SoCavg= socData(timestep);
    OCVcell=interpn(OcvLuts{1}.Dims.soc,OcvLuts{1}.Dims.temp,OcvLuts{1}.Components.ocv,SoCavg,tempData(timestep),'makima');% ocv at avearge soc
    OCVP_surf=interpn(OcvLuts{2}.Dims.soc,OcvLuts{2}.Dims.temp,OcvLuts{2}.Components.U_p,1-SoCPs,tempData(timestep),'makima'); %ocv at surface soc
    OCVN_surf=interpn(OcvLuts{2}.Dims.soc,OcvLuts{2}.Dims.temp,OcvLuts{2}.Components.U_n,SoCNs,tempData(timestep),'makima');
    %Vdiff(timestep)=-kd*(OCVcell-(OCVP_surf-OCVN_surf));
    %v_Out(timestep)=IR0+Vrc+Vdiff(timestep)+OCVcell+U_hyst;
    v_Out(timestep)=IR0+Vrc+U_hyst+kd*(OCVP_surf-OCVN_surf);
    end

    voltOut=v_Out';

    % for i=2:length(timeData)
    %     
    % voltOut(i)=v_Out(N*(i-1)+1);
    % end
    % voltOut(1)=v_Out(1);


end

