% ==============================================================================
% This routine considers the optimal control problem that controls the velocity 
% of the melting interface in cell thawing.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Preparation
% Add paths
addpath(genpath('Saved Data'),'Python','Events','Input Data','Plotting', ...
    'Objective Functions','Optimal Control', 'PDEs','Calculations')
% Add path for CaSADi here.


% Typical input data
input_def = get_inputdata;  % unpack all default input data
input_def.nc = 2;  % control interval
input = input_processing(input_def);  % processing the input  
Tb_def = input.Tb;  % default heater temperature 
n1 = input.n1;
endtime = input.endtime_Itf;  % final time
tspan = [(0:input.dt_plot:endtime)'; endtime];  % time span with interval = dt_plot
tspan = unique(tspan);  % remove duplicated value
T_ini = input.IC;  % initial condition
dSdt_obj = input.dSdt_obj;  % target dSdt
S_target = 1+tspan*dSdt_obj;  % target interface position

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

filename = 'Tb_MeltingSpeed_nc=12_fmincon.mat';
Data_opt = load(filename).Data;

% Solve the ODEs with the optimal profiles
[t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
S_opt = T_opt(:,end);
S_target = 1+t_opt*dSdt_obj;
RMSE = sqrt(sum((S_opt-S_target).^2)/length(t_opt));
e2 = cal_dSdt_2norm(t_opt,S_opt,dSdt_obj);
disp(['e2 = ' , num2str(e2)])

figure; plot_interface(t_opt,S_opt) 
figure; plot_Tb(Data_opt.tb,input.temp_non(Data_opt.Tb_opt))
figure; plot_dSdt(t_opt,S_opt,input)

end


%% Opimization using Ipopt with a piecewise linear control
switch Ipopt
case 'on'
disp('Solving the optimal control problem with Ipopt')

% Solve the optimal control problem
Data = OCP_meltingspeed_Ipopt(input,'');

% Solve the ODEs with the optimal profiles
[t_Ipopt,T_Ipopt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_Ipopt = T_Ipopt(:,end);
Data.RMSE = sqrt(sum((S_Ipopt-S_target).^2)/length(t_Ipopt));
Data.e2 = cal_dSdt_2norm(t_Ipopt,S_Ipopt,dSdt_obj);
disp(['e2_Ipopt = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_Ipopt,S_Ipopt,input)

end


%% Opimization using fmincon with a piecewise linear control
switch fmincon
case 'on'
disp('Solving the optimal control problem with fmincon')

% Solve the optimal control problem
Data = OCP_meltingspeed_fmincon(input,'');

% Solve the ODEs with the optimal profiles
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Data.RMSE = sqrt(sum((S_fmincon-S_target).^2)/length(t_fmincon));
Data.e2 = cal_dSdt_2norm(t_fmincon,S_fmincon,dSdt_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_fmincon,S_fmincon,input)

end


%% Direct single shooting with CasADi
switch CasADi_ss
case 'on'
disp('Solving the optimal control problem with CasADi, single shooting')

% Solve the optimal control problem
Data = OCP_meltingspeed_CSD_ss(input,'');

% Solve the ODEs with the optimal profiles
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
S_target = 1+t_CSD*dSdt_obj;
Data.RMSE = sqrt(sum((S_CSD-S_target).^2)/length(t_CSD));
Data.e2 = cal_dSdt_2norm(t_CSD,S_CSD,input.dSdt_obj);

end


%% Solving the optimal control problem with the DAE-based technique
switch DAE
case 'on'
disp('Solving the optimal control problem with the DAE-based technique')

% Load the Python file
py_file = 'CellThawing_DAE_Interface.py';

% DAE Solver from Python
Tb0 = 0.2186;  % consistent initial condition
%Tb0 = 0.5;  % inconsistent initial condition
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB',x=int64(1),y=dSdt_obj,z=Tb0);
Theta_opt_DAE = double(output_py{2})';
Data.Tb_opt = input.temp_dim(Theta_opt_DAE);
Data.tb = double(output_py{3})';
Data.twall = output_py{1};
Data.ext = '';
Data.Tb_ini = Tb0;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
    
% Solve the ODEs with the optimal profiles from the DAE-based technique
[t_DAE,T_DAE] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data.Tb_opt,Data.tb],input), tspan, T_ini, option_ode2);
S_DAE = T_DAE(:,end);
Data.RMSE = sqrt(sum((S_DAE-S_target).^2)/length(t_DAE));
Data.e2 = cal_dSdt_2norm(t_DAE,S_DAE ,dSdt_obj);
disp(['e2_DAE = ' , num2str(Data.e2)])

figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_DAE,S_DAE,input)

end
