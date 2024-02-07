function outputs = cal_dT_2norm(t,T,Tb,dT_obj)

nt = length(t);
e2 = (Tb-T-dT_obj);
outputs = (1/sqrt(nt))*norm(e2,2);

return