%% Load experiment and OCV lookup 

clear;
load CC_25.mat;
BrOcv = gdParam.OCV_Fill_Sparse_OCV("J:\01_Cell_Database\Cells\Samsung\48X\OCV\HysteresisFull\Rev_1\48X_HysteresisFull_1001z_7T.mat");
s48xFolder = 'J:\01_Cell_Database\Cells\Samsung\48X\ECN\2RC\Rev_1\';
BrEcnName = '48X_2RC_21z_7T_5I_14.mat';
BrEcn = gdParam.ECN_Fill_Sparse_ECN([s48xFolder BrEcnName]);
load lincc_25.mat;

k0=[1000;  0.7]; % tauD, kd, Ea
lowBound = [  500;   k0(2)*0.5];
upBound =  [  1800;  k0(2)*2];

iniPopSpread=5;
initGaPopSize = 70;
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
				initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',1000.*length(k0),...
				'HybridFcn',{@fmincon,hybridopts});


 %% N=5, 15A

N=5;
currData=lincc_25{N}(4:1100,1);
socData=lincc_25{N}(4:1100,3);
voltageData=lincc_25{N}(4:1100,2);
tempData=lincc_25{N}(4:1100,4);
timeData=4:length(currData)+3;
dt=1;

genAlgOpts.MaxTime = 60*5;
maxOptTime = 60*5;

%  lowBound(7:8)= k0(7:8);
%  upBound (7:8)= k0(7:8);

globObj = @(k)ECN_Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,BrOcv,BrEcn);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT_lincc(:,2)=xOptTmp;
temp=xOptTmp;
save param_lincc15.mat xOptTmp
%% N=1, 12.5A
% N=1;
% currData=lincc_25{N}(4:2000,1);
% socData=lincc_25{N}(4:2000,3);
% voltageData=lincc_25{N}(4:2000,2);
% tempData=lincc_25{N}(4:2000,4);
% timeData=4:length(currData)+3;
% dt=1;
% 
% k0=xOptTmp;
% lowBound =[k0(1)*0.8;  k0(2)*0.5;    k0(3)];
% upBound = [k0(1)*1.2;   k0(2)*2;    k0(3)];
% globObj = @(k)ECN_Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,BrOcv,BrEcn);
% [xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);
% 
% 
% LUT_lincc(:,1)=xOptTmp;
% save param_lincc125.mat xOptTmp


%% N=8 20A

% N=8;
% currData=lincc_25{N}(4:600,1);
% socData=lincc_25{N}(4:600,3);
% voltageData=lincc_25{N}(4:600,2);
% tempData=lincc_25{N}(4:600,4);
% timeData=4:length(currData)+3;
% dt=1;
% 
% k0=temp;
% lowBound =[k0(1)*0.5;  k0(2)*0.5;    k0(3)];
% upBound = [k0(1)*2;   k0(2)*2;    k0(3)];
% globObj = @(k)ECN_Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,BrOcv,BrEcn);
% [xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);
% 
% 
% LUT_lincc(:,3)=xOptTmp;
% save param_lincc20.mat xOptTmp



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

error_fit=mean(abs(Vsim-voltageData))
 

%% save lut



% 
 save LUT_lincc.mat LUT_lincc;



 %% Test run
% % 
N=9;
params=xOptTmp;

currData_t=lincc_25{N}(4:1100,1);
socData_t=lincc_25{N}(4:1100,3);
voltageData_t=lincc_25{N}(4:1100,2);
tempData_t=lincc_25{N}(4:1100,4);
timeData=4:length(currData_t)+3;
% 
[Vsim]=ECN_diffusion_model(params,currData_t,timeData,socData_t,tempData_t,BrOcv,BrEcn);
error_val=mean(abs(Vsim-voltageData_t))
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
