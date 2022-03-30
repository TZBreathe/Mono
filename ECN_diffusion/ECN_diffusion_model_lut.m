%% A simple diffusion model

function voltOut=ECN_diffusion_model_lut(k,currData,timeData,socData,tempData,OcvLuts,ECNParams)
%params

N=1; % controls timestep
% dt=timeData(2)-timeData(1);
dt=1;
Tref=298;
% Crate=[0.1; 0.5; 1; 2; 3; 4; 5].*4.75; %CC charge
Crate=[12.5; 15; 20]; % Lincc

SoC0=socData(1);
voltOut=ones(length(timeData),1);
Vrc1=0;Vrc2=0;

% diffusion settings
Q=4.7455*3600; %capacity, should be an argument but lazy for the moment
R = 1; % particle radius [m]
Nr = 10; % number of "shells" radially
dR = R/Nr; % width of each "shell"
Sa = 4*pi*(R*(1:Nr)/Nr).^2; % outer surface area of each shell
dV = (4/3)*pi*((R*(1:Nr)/Nr).^3-(R*(0:Nr-1)/Nr).^3); % vol. of ea. shell
SoC = SoC0*ones(1,(Nr)); % concentration profile versus "r" dimension
SoCs = zeros(size(timeData)); % concentration at surface
SoCavg=SoC0*ones(size(timeData));
SoCs(1) = SoC0;

h(1)=-1;
% h(1)=1;
k_hyst=10;

%if using finer timestep 'interploate' input current array to larger size
times=N*length(timeData);
curr_Data=currData;

    tau_0=interp1(Crate,k(1,:),max(curr_Data),'linear',k(1,3));
    kd=interp1(Crate,k(2,:),max(curr_Data),'linear',k(2,3));
    Ea=interp1(Crate,k(3,:),max(curr_Data),'linear',k(3,3));
% for ii=2:length(currData)
%   curr_Data(N*ii+1:N*ii+10)=currData(ii);
% end
% curr_Data(1:N)=currData(1);

% calc diffusion 
for timestep = 1:times

    
    tau=tau_0*exp(-Ea/8.314*(-1/(273+tempData(timestep))+1/Tref));  
    R0=interpn(ECNParams.Dims.soc, ECNParams.Dims.temp, ECNParams.Dims.crate,ECNParams.Components.R_0,socData(timestep),tempData(timestep),curr_Data(timestep)./4.7,'linear',ECNParams.Components.R_0(ceil(socData(timestep)*21),6,8));
    R1=interpn(ECNParams.Dims.soc, ECNParams.Dims.temp, ECNParams.Dims.crate,ECNParams.Components.R_1,socData(timestep),tempData(timestep),curr_Data(timestep)./4.7,'linear',ECNParams.Components.R_1(ceil(socData(timestep)*21),6,8));
    R2=interpn(ECNParams.Dims.soc, ECNParams.Dims.temp, ECNParams.Dims.crate,ECNParams.Components.R_2,socData(timestep),tempData(timestep),curr_Data(timestep)./4.7,'linear',ECNParams.Components.R_2(ceil(socData(timestep)*21),6,8));
    C1=interpn(ECNParams.Dims.soc, ECNParams.Dims.temp, ECNParams.Dims.crate,ECNParams.Components.C_1,socData(timestep),tempData(timestep),curr_Data(timestep)./4.7,'linear',ECNParams.Components.C_1(ceil(socData(timestep)*21),6,8));
    C2=interpn(ECNParams.Dims.soc, ECNParams.Dims.temp, ECNParams.Dims.crate,ECNParams.Components.C_2,socData(timestep),tempData(timestep),curr_Data(timestep)./4.7,'linear',ECNParams.Components.C_2(ceil(socData(timestep)*21),6,8));
    
    IR0=R0.*curr_Data(timestep); 
    Vrc1=R1*(exp(-dt/R1/C1).*(Vrc1/R1)+(1-exp(-dt/R1/C1)).*curr_Data(timestep));
    Vrc2=R2*(exp(-dt/R2/C2).*(Vrc2/R2)+(1-exp(-dt/R2/C2)).*curr_Data(timestep));

    M_hyst=interpn(OcvLuts.Dims.soc, OcvLuts.Dims.temp,OcvLuts.Components.hystAmp,socData(timestep),tempData(timestep),'makima'); 
    M0=interpn(OcvLuts.Dims.soc, OcvLuts.Dims.temp,OcvLuts.Components.hystInst,socData(timestep),tempData(timestep),'makima');
    h=exp(-dt*abs(curr_Data(timestep))*k_hyst/(Q))*h+sign(curr_Data(timestep))*(1-exp(-dt*abs(curr_Data(timestep)*k_hyst/(Q))));
    U_hyst=M_hyst.*h+sign(curr_Data(timestep)).*M0;

    flux = -1/tau*diff(SoC)/dR; % flux at surfaces between "bins"
    M = flux.*Sa(1:end-1); % total SoC crossing surface between bins
    SoC= SoC+ ([0 M] - [M 0])*dt./dV; % conc. change via diffusion
    SoC(end) = SoC(end) + (curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end); % at boundary
    SoCs = min(1,SoC(end)); % surface soc
    SoCavg= socData(timestep);
    OCVcell=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.ocv,SoCavg,tempData(timestep),'linear',OcvLuts.Components.ocv(ceil(SoCavg*1000),7));% ocv at avearge soc
    %OCVcell=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.ocv,SoCavg(timestep),tempData(timestep),'makima');
    OCVcell_surf=interpn(OcvLuts.Dims.soc,OcvLuts.Dims.temp,OcvLuts.Components.ocv,SoCs,tempData(timestep),'linear',OcvLuts.Components.ocv(max(ceil(SoCs*1000),7),1)); %ocv at surface soc
    Vdiff=-kd*(OCVcell-OCVcell_surf);
    v_Out(timestep)=IR0+Vrc1+Vrc2+Vdiff+OCVcell+U_hyst;
%    v_Out(timestep)=IR0+Vrc1+Vrc2+OCVcell+U_hyst;

end

    voltOut=v_Out';

end

