function outputs = obj_inputpower_pl(Tb, tb, tspan, dT_obj, T_ini, option, input)

    Heater_Profile = [Tb, tb];
    [t,T] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Heater_Profile,input), tspan, T_ini, option);
   
    Tavg = cal_Tavg(T(:,input.n1+1:end-1));
    Tb_tg = input.temp_non(interp1(tb, Tb, t));
    outputs = sum((Tb_tg-Tavg-dT_obj).^2);

end