%% Load experiment and OCV lookup 

clear;
BrOcv = gdParam.OCV_Fill_Sparse_OCV("J:\01_Cell_Database\Cells\Samsung\48X\OCV\HysteresisFull\Rev_1\48X_HysteresisFull_1001z_7T.mat");
s48xFolder = 'J:\01_Cell_Database\Cells\Samsung\48X\ECN\2RC\Rev_1\';
BrEcnName = '48X_2RC_21z_7T_9I.mat';
% BrEcn = gdParam.ECN_Fill_Sparse_ECN([s48xFolder BrEcnName]);
 load([s48xFolder BrEcnName]);
load lincc_35.mat;

k0=[1000;  0.7; 8000]; % tauD, kd, Ea
lowBound = [500;   k0(2)*0.5; 3000];
upBound =  [2800;  k0(2)*2; 15000];

iniPopSpread=7;
initGaPopSize = 50;
randInitPopFactor = max(min(randn(length(k0),initGaPopSize)/2,1),-1);
initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
genAlgOpts.InitialPopulationMatrix = initGaPop';		
genAlgOpts.MaxTime = 60*5;

maxOptTime = 60*5; % Set optimisation time
maxHybridFminIters = 300; % Max number of IPA iterations


% IPA optimisation options
hybridopts = optimoptions('fmincon','Display','none','MaxIterations',maxHybridFminIters,'UseParallel',true);
				
				% GA optimisation options
genAlgOpts = optimoptions('ga','Display','iter','InitialPopulationMatrix',...
				initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',2000.*length(k0),...
				'HybridFcn',{@fmincon,hybridopts});


 %% N=3, 20A

N=3;
currData=lincc_35{N}(4:1000,1);
socData=lincc_35{N}(4:1000,3);
voltageData=lincc_35{N}(4:1000,2);
tempData=lincc_35{N}(4:1000,4);
timeData=4:length(currData)+3;
dt=1;

genAlgOpts.MaxTime = 60*5;
maxOptTime = 60*5;

%  lowBound(7:8)= k0(7:8);
%  upBound (7:8)= k0(7:8);

globObj = @(k)ECN_Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,BrOcv,BrEcn);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT_lincc_35(:,2)=xOptTmp;
temp=xOptTmp;
% save param_lincc15.mat xOptTmp
%% N=1, 12.5A
N=1;
currData=lincc_35{N}(4:1100,1);
socData=lincc_35{N}(4:1100,3);
voltageData=lincc_35{N}(4:1100,2);
tempData=lincc_35{N}(4:1100,4);
timeData=4:length(currData)+3;
dt=1;

k0=xOptTmp;
lowBound =[k0(1)*0.8;  k0(2)*0.5;    k0(3)*0.7];
upBound = [k0(1)*1.2;   k0(2)*2;    k0(3)*1.4];
globObj = @(k)ECN_Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,BrOcv,BrEcn);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);


LUT_lincc_35(:,1)=xOptTmp;
% save param_lincc125.mat xOptTmp




%% Run and plot
params=xOptTmp;

[Vsim]=ECN_diffusion_model(params,currData,timeData,socData,tempData,BrOcv,BrEcn);

figure();
hold on;
plot(socData,Vsim,'bl');
plot(socData,voltageData);
xlabel('SoC');
ylabel('Voltage');

yyaxis right
plot(socData,currData,'red');
ylabel('Current','color','red');
legend('Model','Exp','Current','location','southeast');
hold off;

% error_fit=mean(abs(Vsim-voltageData))
error_fit=sqrt(mean((Vsim-voltageData).^2))
 

%% save lut
 save LUT_lincc_35.mat LUT_lincc_35;

 %% Test run

N=2;
% params=xOptTmp;
load('LUT_lincc_35.mat');
params=LUT_lincc_35;

currData_t=lincc_35{N}(4:1000,1);
socData_t=lincc_35{N}(4:1000,3);
voltageData_t=lincc_35{N}(4:1000,2);
tempData_t=lincc_35{N}(4:1000,4);
timeData=4:length(currData_t)+3;
% 
[Vsim]=ECN_diffusion_model_lut(params,currData_t,timeData,socData_t,tempData_t,BrOcv,BrEcn);
 error_val=sqrt(mean((Vsim-voltageData_t).^2))
% mean(abs(Vsim-voltageData_t))
% 
% 
% 
figure();
hold on;
plot(socData_t,Vsim,'bl');
plot(socData_t,voltageData_t);
xlabel('SoC');
ylabel('Voltage');

yyaxis right
ax = gca;
ax.YColor = 'r';
plot(socData_t,currData_t,'red');
ylabel('Current','color','red');
legend('Model','Exp','Current','location','southeast');
hold off;



% plot(socData_t,voltageData_t-Vsim);
% 
% % r=linspace(0,1,20);
% % hold on;
% % for i=1:200:1000
% %     plot(r,socr(i,:));
% % end
% % xlabel('r');
% % ylabel('SOC');
% % legend('1s', '200s','400s','600s','8000s');
% % hold off;
