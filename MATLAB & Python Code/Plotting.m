% ==============================================================================
% This is a routine for plotting.
%
% Created by Prakitr Srisuma, 
% PhD in Computational Science and Engineering, MIT
% ==============================================================================
close all; clear; clc;

%% Preparation
% Add paths
addpath(genpath('Saved Data'),'Python','Events','Input Data','Plotting', ...
    'Objective Functions','Optimal Control', 'PDEs','Calculations','Exporting Graphics','Figures')
addpath(['C:\Users\pkta1\OneDrive - Chulalongkorn University\Academic\MIT\PhD_CSE\Research\' ...
    'Cell Thawing\Optimal Control\MATLAB\CasADi\casadi-3.6.3-windows64-matlab2018b'])

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

% ODE solvers
tol_ode = input.tol_ode;
option_ode1 = input.option_ode1;
option_ode2 = input.option_ode2;

% Optimal control
nrun = input.nrun;  % number of simulation runs
nc = input.nc;  % number of control intervals

% Option
Prob3_revised = 'on';
Prob5 = 'on';


%% Plotting Problem 2 ACC Paper
switch Prob3_revised
case 'on'

file = {'Tb_MeltingSpeed_nc=12_Ipopt.mat','Tb_MeltingSpeed_nc=12_fmincon.mat',...
    'Tb_MeltingSpeed_nc=12_CSD_ss.mat','Tb_MeltingSpeed_DAE'};
method = {'opt\_Ipopt','opt\_fmincon','opt\_CasADi', 'sim\_DAE'};
linspec = {'-x','-^','-square','-o','-diamond','-*'};
msize = [7;5;6;5];
lw = [2;1;1;1];
color = {'b',[0 0.62 0.17],[0.6350 0.0780 0.1840],[0.9290 0.5940 0.1250]};

N = 4;
nplot = N*4*12;
tplot = linspace(0,7,nplot)';

figure
fig = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
for i = 1:4
    nexttile(1)
    Data = load(file{i}).Data;
    Tb = input.temp_non(Data.Tb_opt);
    tb = Data.tb;
    yplot = interp1(tb,Tb,tplot);

    if i == 1|| i == 2
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('\Theta_{\itb}')
    xlabel('\tau')
    xlim([0 7])
    xticks(0:1:7)
    legend('location','southeast')
    hold on

end
    graphics_setup('1by3.5_2')
    grid on
    text(.02,0.95,'(A)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 
    nexttile(2)

for i = 1:4
    Data_opt = load(file{i}).Data;
    [t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
    S_opt = T_opt(:,end);
    dSdt = cal_dSdt(t_opt,S_opt);
    dSdt = [dSdt; dSdt(end)];
    tb = t_opt;
    yplot = interp1(tb,dSdt,tplot);

    if i == 1|| i == 2
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[1;((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('{\itdS}/{\itd}\tau')
    xlabel('\tau')
    xlim([0 7])
    ylim([-0.104 -0.096])
    xticks(0:1:7)
    legend('location','southeast')
    hold on
end
    graphics_setup('1by3.5_2')
    grid on
    text(.02,0.95,'(B)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 
    nexttile(3)

    
for i = 1:4
    Data_opt = load(file{i}).Data;
    [t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
    S_opt = T_opt(:,end);
    tb = t_opt;
    yplot = interp1(tb,S_opt,tplot);

    if i == 1
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('{\itS}')
    xlabel('\tau')
    xlim([0 7])
    xticks(0:1:7)
    legend('location','best')
    hold on
end
    graphics_setup('1by3.5_2')
    grid on
    text(.15,0.95,'(C)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 

end


%% Plotting Problem 1 ACC Paper
switch Prob5
case 'on'

file = {'Tb_InputPower_nc=12_Ipopt.mat','Tb_InputPower_nc=12_fmincon.mat',...
    'Tb_InputPower_nc=12_CSD_ss.mat','Tb_InputPower_DAE'};
method = {'opt\_Ipopt','opt\_fmincon','opt\_CasADi', 'sim\_DAE'};
linspec = {'-x','-^','-square','-o','-diamond','-*'};
msize = [7;5;6;5];
lw = [2;1;1;1];
color = {'b',[0 0.62 0.17],[0.6350 0.0780 0.1840],[0.9290 0.5940 0.1250]};

N = 4;
nplot = N*4*12;
tplot = linspace(0,7,nplot)';

figure
fig = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
for i = 1:4
    nexttile(1)
    Data = load(file{i}).Data;
    Tb = input.temp_non(Data.Tb_opt);
    tb = Data.tb;
    yplot = interp1(tb,Tb,tplot);

    if i == 1
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('\Theta_{\itb}')
    xlabel('\tau')
    xlim([0 7])
    xticks(0:1:7)
    legend('location','southeast')
    hold on

end
    graphics_setup('1by3.5')
    grid on
    text(.02,0.95,'(A)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 
    nexttile(2)

    
for i = 1:4
    Data_opt = load(file{i}).Data;
    [t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
    S_opt = T_opt(:,end);
    Tavg = cal_Tavg(T_opt(:,n1+1:end-1));
    tb = t_opt;
    yplot = interp1(tb,Tavg,tplot);

    if i == 1
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('{\Theta_{avg}}')
    xlabel('\tau')
    xlim([0 7])
    xticks(0:1:7)
    legend('location','best')
    hold on
end
    graphics_setup('1by3.5')
    grid on
    text(.02,0.95,'(B)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 
    nexttile(3)

for i = 1:4
    Data_opt = load(file{i}).Data;
    [t_opt,T_opt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,[Data_opt.Tb_opt,Data_opt.tb],input), tspan, T_ini, option_ode2);
    Tavg = cal_Tavg(T_opt(:,n1+1:end-1));
    tb = t_opt;
    Tavg = interp1(tb,Tavg,tplot);
    Tb = input.temp_non(Data.Tb_opt);
    Tb = interp1(tb, Tb, tplot);
    yplot = Tb - Tavg;
    
    if i == 1
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',[((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))';length(tplot)]);
    else
        plot(tplot,yplot,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'DisplayName',method {i},'linewidth',lw(i),'markersize',msize(i),'MarkerIndices',((1 + ((nplot/N)/4)*(i-1)):nplot/N:length(tplot))');
    end

    ylabel('\Theta_{\itb}–\Theta_{avg}')
    xlabel('\tau')
    xlim([0 7])
    ylim([0 0.2])
    xticks(0:1:7)
    legend('location','southeast')
    hold on
end
    graphics_setup('1by3.5_2')
    grid on
    text(.02,0.95,'(C)','Units','normalized','FontSize', 12 ,'fontweight', 'bold') 

end
