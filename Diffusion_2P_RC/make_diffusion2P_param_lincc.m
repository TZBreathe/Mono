%% Load experiment and OCV lookup 

clear;
load OCV;
load lincc_25.mat;
ocvData=OCV;
load LUT25.mat


k0=[0.02; 0.01; 100/0.01; 800; 2000; 1.0; 12000; 9000]; %R0, R1,  C1, tauP, tauD,, kd, Ea1, Ea2
lowBound = [k0(1)/3;  k0(2)/3;  k0(3)/3;   500;  800;  k0(6)*0.8;  6000;  6000];
upBound =  [k0(1)*3;  k0(2)*3;  k0(3)*3;  2500;  3000; k0(6)*1.2;  18000; 14000];

iniPopSpread=4;
initGaPopSize = 70;
randInitPopFactor = max(min(randn(length(k0),initGaPopSize)/2,1),-1);
initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
genAlgOpts.InitialPopulationMatrix = initGaPop';		
genAlgOpts.MaxTime = 60*10;

maxOptTime = 60*10; % Set optimisation time
maxHybridFminIters = 300; % Max number of IPA iterations


% IPA optimisation options
hybridopts = optimoptions('fmincon','Display','none','MaxIterations',maxHybridFminIters,'UseParallel',true);
				
				% GA optimisation options
genAlgOpts = optimoptions('ga','Display','iter','InitialPopulationMatrix',...
				initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',1000.*length(k0),...
				'HybridFcn',{@fmincon,hybridopts});
ocvData=OCV;

 %% N=5, 15A

N=5;
currData=lincc_25{N}(4:1100,1);
socData=lincc_25{N}(4:1100,3);
voltageData=lincc_25{N}(4:1100,2);
tempData=lincc_25{N}(4:1100,4);
timeData=4:length(currData)+3;
dt=1;

genAlgOpts.MaxTime = 60*15;
maxOptTime = 60*15;

%  lowBound(7:8)= k0(7:8);
%  upBound (7:8)= k0(7:8);

globObj = @(k)Diffusion2P_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData,dt);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT_lincc(:,2)=xOptTmp;
temp=xOptTmp;
save param_lincc15.mat xOptTmp
%% N=1, 12.5A
N=1;
currData=lincc_25{N}(4:2000,1);
socData=lincc_25{N}(4:2000,3);
voltageData=lincc_25{N}(4:2000,2);
tempData=lincc_25{N}(4:2000,4);
timeData=4:length(currData)+3;
dt=1;

k0=xOptTmp;
lowBound = [k0(1)*0.8;  k0(2)*0.7;    k0(3)*0.7;   k0(4)*0.8;  k0(5)*0.8;  k0(6)*0.8; k0(7); k0(8)];
upBound = [k0(1)*1.2;   k0(2)*1.4;    k0(3)*1.4;   k0(4)*1.2;  k0(5)*1.2;  k0(6)*1.2; k0(7); k0(8)];
globObj = @(k)Diffusion2P_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData,dt);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);


LUT_lincc(:,1)=xOptTmp;
save param_lincc125.mat xOptTmp


%% N=8 20A

N=8;
currData=lincc_25{N}(4:600,1);
socData=lincc_25{N}(4:600,3);
voltageData=lincc_25{N}(4:600,2);
tempData=lincc_25{N}(4:600,4);
timeData=4:length(currData)+3;
dt=1;

k0=temp;
lowBound = [k0(1)*0.8;  k0(2)*0.7;    k0(3)*0.7;   k0(4)*0.7;  k0(5)*0.7;  k0(6)*0.8; k0(7); k0(8)];
upBound = [k0(1)*1.2;   k0(2)*1.4;    k0(3)*1.4;   k0(4)*1.4;  k0(5)*1.4;  k0(6)*1.2; k0(7); k0(8)];
globObj = @(k)Diffusion2P_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData,dt);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);


LUT_lincc(:,3)=xOptTmp;
save param_lincc20.mat xOptTmp



%% Run and plot
params=xOptTmp;

[Vsim]=diffusion2P_model(params,currData,timeData,socData,tempData,ocvData,1);

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
N=6;
load LUT_lincc.mat
load lincc_25.mat;
load OCV;
ocvData=OCV;
params=LUT_lincc;
currData_t=lincc_25{N}(4:1100,1);
socData_t=lincc_25{N}(4:1100,3);
voltageData_t=lincc_25{N}(4:1100,2);
tempData_t=lincc_25{N}(4:1100,4);
timeData=4:length(currData_t)+3;
% 
Vsim=diffusion2P_model_run_lut(params,currData_t,timeData,socData_t,tempData_t,ocvData);
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
