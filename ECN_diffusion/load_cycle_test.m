%%

clear;
BrOcv = gdParam.OCV_Fill_Sparse_OCV("J:\01_Cell_Database\Cells\Samsung\48X\OCV\HysteresisFull\Rev_1\48X_HysteresisFull_1001z_7T.mat");
s48xFolder = 'J:\01_Cell_Database\Cells\Samsung\48X\ECN\2RC\Rev_1\';
BrEcnName = '48X_2RC_21z_7T_9I.mat';
% BrEcn = gdParam.ECN_Fill_Sparse_ECN([s48xFolder BrEcnName]);
 load([s48xFolder BrEcnName]);

%% normal disch profile
load('C:\Users\Teng\Breathe Battery Technologies\Gold Dust - Documents\03_Load_Cycles\Cycles\Car\Supercar\Discharge\Teide_normal_Supercar_Discharge_1_20210412_14-16-35.mat');

time_len= RunData.dataTable{end-1,'timeTest'}-RunData.dataTable{1,'timeTest'};
timeData=linspace(0,time_len,time_len/1);

currData=RunData.dataTable{1:end,'currCell'};
currData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'currCell'},timeData);
voltageData=RunData.dataTable{1:1:end,'voltCell'};
voltData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'voltCell'},timeData);

% soc0=interp1(BrOcv.Components.ocv(:,5),BrOcv.Dims.soc,voltageData(1));
soc0=1;
ahData=RunData.dataTable{1:1:end,'ahTotal'};
ahData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'ahTotal'},timeData);
socData_intp=soc0+ahData_intp/4.75;

tempData=RunData.dataTable{1:1:end,'tempCell'};
tempData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'tempCell'},timeData);

%  hold on
% plot(timeData,ahData_intp)
% plot(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'ahTotal'},'o');
%%
load('LUT_lincc');
params=LUT_lincc(:,1);

[Vsim,h]=ECN_diffusion_model_LC(params,currData_intp,timeData,socData_intp,tempData_intp,BrOcv,BrEcn);

figure();
hold on;
plot(timeData,Vsim,'bl');
plot(timeData,voltData_intp);
xlabel('SoC');
ylabel('Voltage');

RMSE=sqrt(mean((voltData_intp'-Vsim).^2));
Err_max=max(abs((voltData_intp'-Vsim)));

%% demanding
load('C:\Users\Teng\Breathe Battery Technologies\Gold Dust - Documents\03_Load_Cycles\Cycles\Car\Supercar\Discharge\Teide_demanding_Supercar_Discharge_1_20210412_14-16-32.mat');

time_len= RunData.dataTable{end-1,'timeTest'}-RunData.dataTable{1,'timeTest'};
timeData=linspace(0,time_len,time_len/0.2);

currData=RunData.dataTable{1:1:end,'currCell'};
currData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'currCell'},timeData);
voltageData=RunData.dataTable{1:1:end,'voltCell'};
voltData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'voltCell'},timeData);

soc0=interp1(BrOcv.Components.ocv(:,5),BrOcv.Dims.soc,voltageData(1));
ahData=RunData.dataTable{1:1:end,'ahTotal'};
ahData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'ahTotal'},timeData);
socData_intp=soc0+ahData_intp/4.75;

tempData=RunData.dataTable{1:1:end,'tempCell'};
tempData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'tempCell'},timeData);
% 
% hold on
% plot(timeData,ahData_intp)
% plot(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'ahTotal'},'o');

%%
load('LUT_lincc');
params=LUT_lincc(:,3);
[Vsim,socs]=ECN_diffusion_model_LC(params,currData_intp,timeData,socData_intp,tempData_intp,BrOcv,BrEcn);

figure();
hold on;
plot(timeData,Vsim,'bl');
plot(timeData,voltData_intp);
xlabel('SoC');
ylabel('Voltage');

RMSE=sqrt(mean((voltData_intp'-Vsim).^2));
Err_max=max(abs((voltData_intp'-Vsim)));

%% very demanding

load('C:\Users\Teng\Breathe Battery Technologies\Gold Dust - Documents\03_Load_Cycles\Cycles\Car\Supercar\Discharge\Teide_very_demanding_Supercar_Discharge_1_20210412_14-16-23.mat');

time_len= RunData.dataTable{end-1,'timeTest'}-RunData.dataTable{1,'timeTest'};
timeData=linspace(0,time_len,time_len/0.1);

currData=RunData.dataTable{1:1:end,'currCell'};
currData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'currCell'},timeData);
voltageData=RunData.dataTable{1:1:end,'voltCell'};
voltData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'voltCell'},timeData);

soc0=interp1(BrOcv.Components.ocv(:,5),BrOcv.Dims.soc,voltageData(1));
ahData=RunData.dataTable{1:1:end,'ahTotal'};
ahData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'ahTotal'},timeData);
socData_intp=soc0+ahData_intp/4.75;

tempData=RunData.dataTable{1:1:end,'tempCell'};
tempData_intp=interp1(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'tempCell'},timeData);
% 
hold on
plot(timeData,voltData_intp)
plot(RunData.dataTable{:,'timeTest'}-RunData.dataTable{1,'timeTest'},RunData.dataTable{1:1:end,'voltCell'});

%%
load('LUT_lincc');
params=LUT_lincc(:,3);
[Vsim,socs]=ECN_diffusion_model_LC(params,currData_intp,timeData,socData_intp,tempData_intp,BrOcv,BrEcn);

figure();
hold on;
plot(timeData,Vsim,'bl');
plot(timeData,voltData_intp);
xlabel('SoC');
ylabel('Voltage');

RMSE=sqrt(mean((voltData_intp'-Vsim).^2));
Err_max=max(abs((voltData_intp'-Vsim)));