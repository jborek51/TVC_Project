%%  Obtain Motor Properties 
num = xlsread('Motor_Data.xlsx',1,'B3:D30');
Motor.F_avg = mean(num(:,2));                         %   N - Average thrust 
Motor.F_r = num(2,3);                                 %   N - Rail thrust
Motor.F_th = num(:,2);                                %   N - Thrust vector 
Motor.I = num(1,3);                                   %   Ns - Total impulse
Motor.m_f = num(4,3);                                 %   kg - Fuel mass
Motor.m_m = num(3,3);                                 %   kg - Motor mass
Motor.t_b = num(5,3);                                 %   s - Burn time
Motor.T_b = num(:,1);                                 %   s - Time vector

%%
save('F15.mat','Motor')
