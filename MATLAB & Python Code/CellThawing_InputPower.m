% ==============================================================================
% This routine considers the optimal control problem that controls the input 
% power (temperature difference) from the heater in cell thawing.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Preparation
% Add paths
addpath(genpath('Saved Data'),'Python','Events','Input Data','Plotting', ...
    'Objective Functions','Optimal Control', 'PDEs','Calculations')
% Add path for CasADi here.

% Typical input data
input_def = get_inputdata;  % unpack all default input data
input_def.nc = 1;  % control interval
input = input_processing(input_def);  % processing the input  
Tb_def = input.Tb;  % default heater temperature 
n1 = input.n1;
endtime = input.endtime_temp;  % final time
tspan = [(0:input.dt_plot:endtime)'; endtime];  % time span with interval = dt_plot
tspan = unique(tspan);  % remove duplicated value
T_ini = input.IC;  % initial condition
dT_obj = input.dT_obj;  % target dT

% ODE solvers
tol_ode = input.tol_ode;
option_ode1 = input.option_ode1;
option_ode2 = input.option_ode2;

% Optimal control
nrun = input.nrun;  % number of simulation runs
nc = input.nc;  % number of control intervals


%% Option selection ('on' or 'off')
orgsim = 'off';  % original simulation
optsim = 'off';  % simulation with the optimal heater temperature
Ipopt = 'off';  % Ipopt
fmincon = 'off';  % fmincon
CasADi_ss = 'off';  % CasADi single shooting
DAE = 'on';  % DAE-based method


%% ODE solver
switch orgsim
case 'on'
[t_sim,T_sim] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,0,input), tspan, T_ini, option_ode2);
S_sim = T_sim(:,end);  % interface position
dTdt = cal_dTdt(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)));

figure; plot_interface(t_sim,S_sim)
figure; plot_Tavg(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)))

end

%% ODE solver
switch optsim
case 'on'

filename = 'Tb_InputPower_nc=12_Ipopt.mat';
Data_opt = load(filename).Data;

% Solve the ODEs with the optimal profiles
[t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
S_opt = T_opt(:,end);
Tavg = cal_Tavg(T_opt(:,input.n1+1:end-1));
Tb = input.temp_non(Data_opt.Tb_opt);
Tb = interp1(Data_opt.tb, Tb, t_opt);
Data.RMSE = sum((Tb-Tavg-dT_obj).^2);
Data.e2 = cal_dT_2norm(t_opt,Tavg,Tb,dT_obj);
disp(['e2_opt = ' , num2str(Data_opt.e2)])

% Plot the result
figure; plot_interface(t_opt,S_opt) 
figure; plot_Tb(Data_opt.tb,input.temp_non(Data_opt.Tb_opt))
figure; plot_dT(t_opt,Tavg,Tb,input)

end


%% Opimization using Ipopt with a piecewise linear control
switch Ipopt
case 'on'
disp('Solving the optimal control problem with Ipopt')

% Solve the optimal control problem
Data = OCP_inputpower_Ipopt(input,'');

% Solve the ODEs with the optimal profiles 
[t_Ipopt,T_Ipopt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_Ipopt = T_Ipopt(:,end);
Tavg = cal_Tavg(T_Ipopt(:,input.n1+1:end-1));
Tb = input.temp_non(Data.Tb_opt);
Tb = interp1(Data.tb, Tb, t_Ipopt);
Data.RMSE = sum((Tb-Tavg-dT_obj).^2);
Data.e2 = cal_dT_2norm(t_Ipopt,Tavg,Tb,dT_obj);
disp(['e2_Ipopt = ' , num2str(Data.e2)])

% Save the result
switch savedata
case 'on'
save(['Saved Data\Tb_InputPower_nc=',num2str(nc),'_Ipopt.mat'],'Data')
end

% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dT(t_Ipopt,Tavg,Tb,input)

end


%% Opimization using fmincon with a piecewise linear control
switch fmincon
case 'on'
disp('Solving the optimal control problem with fmincon')

% Solve the optimal control problem
Data = OCP_inputpower_fmincon(input,'');

% Solve the ODEs with the optimal profiles 
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Tavg = cal_Tavg(T_fmincon(:,input.n1+1:end-1));
Tb = input.temp_non(Data.Tb_opt);
Tb = interp1(Data.tb, Tb, t_fmincon);
Data.RMSE = sum((Tb-Tavg-dT_obj).^2);
Data.e2 = cal_dT_2norm(t_fmincon,Tavg,Tb,dT_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Save the result
switch savedata
case 'on'
save(['Saved Data\Tb_InputPower_nc=',num2str(nc),'_fmincon.mat'],'Data')
end

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dT(t_fmincon,Tavg,Tb,input)

end


%% Direct single shooting with CasADi
switch CasADi_ss
case 'on'
disp('Solving the optimal control problem with CasADi, single shooting')

% Solve the optimal control problem
Data = OCP_inputpower_CSD_ss(input,'');

% Solve the ODEs with the optimal profiles 
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
Tavg = cal_Tavg(T_CSD(:,input.n1+1:end-1));
Tb = input.temp_non(Data.Tb_opt);
Tb = interp1(Data.tb, Tb, t_CSD);
Data.RMSE = sum((Tb-Tavg-dT_obj).^2);
Data.e2 = cal_dT_2norm(t_CSD,Tavg,Tb,dT_obj);
disp(['e2_CSD = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_CSD,S_CSD) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dT(t_CSD,Tavg,Tb,input)

end


%% Solving the optimal control problem with the DAE-based technique
switch DAE
case 'on'
disp('Solving the optimal control problem with the DAE-based technique')

% Run the DAE solver (index 1)
T_ini = [T_ini; input.dT_obj];
M = eye(2*n1+2);
M(end,:) = 0;
options_ode3 = odeset('Events', @(t,T)Event_CompleteMelting(t,T,input),'RelTol',input.tol_ode,'AbsTol',input.tol_ode,'Mass',M);
tic; [t_DAE,T_DAE] = ode15s(@(t,T)DAE_InputPower(t,T,input), tspan, T_ini, options_ode3); 
twall = toc;
S_DAE = T_DAE(:,end-1);  % interface position
Tb = T_DAE(:,end);
Tavg = cal_Tavg(T_DAE(:,input.n1+1:end-2));
Data.Tb_opt = input.temp_dim(Tb);
Data.tb = t_DAE;
Data.twall = twall;
Data.ext = '';
Data.Tb_ini = input.dT_obj;
Data.RMSE = sum((Tb-Tavg-dT_obj).^2);
Data.e2 = cal_dT_2norm(t_DAE,Tavg,Tb,dT_obj);
disp(['e2_DAE = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dT(t_DAE,Tavg,Tb,input)

end
