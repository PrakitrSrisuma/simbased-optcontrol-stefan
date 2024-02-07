function outputs = obj_inputpower_CasADi(t,T,Tb,dT_obj,input)

import casadi.*

Tb_tg = MX.zeros(length(T),1);

for i = 1:length(T)
    dt = t(end)-t(1);
    Tb_tg(i) = (t(i)-t(1))*Tb(2)/(dt) + (t(end)-t(i))*Tb(1)/(dt);
end
Tb_tg = input.temp_non(Tb_tg);

outputs = sum((Tb_tg-T'-dT_obj).^2);


end