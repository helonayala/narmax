
u = load('ts22.dat')'
y = load('ts23.dat')'


ind_est = 1:100;
ind_val = 101:296;

y_est = y(ind_est);
y_val = y(ind_val);

u_est = u(ind_est);
u_val = u(ind_val);
