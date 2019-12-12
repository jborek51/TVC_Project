function [c_ineq,c_eq] = con_fun(U_0,x_0,v_0,t_mode,x_lead,road,...
    m_eff,v_eng_max,veh_param)
%	MPC constraint function used in fmincon

g = veh_param(1);                                   %   m/s^2 - Accel due to gravity
rho = veh_param(2);                                 %   kg/m^3 - Air density
m = veh_param(3);                                   %   kg - Base vehicle mass
Cd = veh_param(4);                                  %   Drag Coefficient
A_ref = veh_param(5);                               %   m^2 - Reference area
Cr = veh_param(6);                                  %   Rolling resistance coefficient
k_D = .5*Cd*rho*A_ref;                              %   kg/m - concentrated drag coef without plattoning

X_vec = road(1,:);                                  %   Road position vector
theta_vec = road(2,:);                              %   Road grade vector
V_lim = road(3,:);                                  %   Road speed limit vector

if t_mode == 1
    t_pred = 24;                                    %   m - Prediction horizon length
    t_ctrl = 16;                                    %   s - Control horizon length
    dt = 1;                                         %   m - Position step
else
    t_pred = 120;                                   %   m - Prediction horizon length
    t_ctrl = 115;                                   %   s - Control horizon length
    dt = 5;                                         %   m - Position step
end
np_steps = t_pred/dt;                               %   # of pred steps
nc_steps = t_ctrl/dt;                               %   # of control steps

X = zeros(1,np_steps);                              %   m - Initialize position vector
s_n = zeros(1,np_steps);                            %   m - Initialize safe following distance vector
ds = zeros(1,np_steps);                             %   m - Initialize following distance vector
V = zeros(1,np_steps);                              %   m/s - Initialize velocity vector
U_app = zeros(1,np_steps);                          %   N - Initialize Ctrl force vector
U_app(1:nc_steps) = U_0;                            %   N - Ctrl force
U_app(nc_steps+1:np_steps) = U_app(nc_steps)...     %   N - Ctrl force
    .*ones(1,np_steps-nc_steps);
theta = zeros(1,np_steps);                          %   rad - Initialize road grade vector
F_grav = zeros(1,np_steps);                         %   N - Initialize gravitational force vector
F_road = zeros(1,np_steps);                         %   N - Initialize rolling resistance force vector
F_aero = zeros(1,np_steps);                         %   N - Initialize aerodynamic force vector          
v_lim = zeros(1,np_steps);                          %   m/s - Initialize max velocity vector

t_n = 2;                                            %   s - Safe following time
r = 1;                                              %   m/s - Rate at which following distance increases
X(1) = x_0;                                         %   m - Initial position
ds(1) = x_lead(1)-x_0;                              %   m - Initial following distance
s_nh = v_0*t_n;                                     %   m - Minimum folowing distance
s_n(1) = min(ds(1),s_nh);                           %   m - Initial safe following distance
V(1) = v_0;                                         %   m/s - Initial velocity

theta(1) = interpone(X_vec,theta_vec,X(1),...         %   rad - Initial road angle
    'linear',theta_vec(end));
v_lim(1) = min(interpone(X_vec,V_lim,X(1),...         %   m/s - Initial max velocity
    'linear',V_lim(end)),v_eng_max);
F_grav(1) = m*g*sin(theta(1));                      %   N - Initial gravitational force
F_road(1) = m*g*Cr*cos(theta(1));                   %   N - Initial rolling resistance force
F_aero(1) = k_D*V(1)^2;                             %   N - Initial aerodynamic force

for i = 2:np_steps    
    V(i) = V(i-1)-dt/m_eff*F_aero(i-1) ...          %   m/s - Next velocity
        -dt/m_eff*F_grav(i-1) ...
        -dt/m_eff*F_road(i-1) ...
        +dt/m_eff*U_app(i-1);
    v_avg = (V(i) + V(i-1))/2;                      %   m/s - Avg velocity over time step
    X(i) = X(i-1) + v_avg*dt;                       %   m - Next position
    ds(i) = x_lead(i)-X(i);                         %   m - Next following distance
    s_n(i) = min(ds(1)+r*(i-1)*dt,s_nh);            %   m - Next safe following distance
    theta(i) = interpone(X_vec,theta_vec,X(i),...     %   rad - Next road angle
        'linear',theta_vec(end));
    v_lim(i) = min(interpone(X_vec,V_lim,X(i),...     %   m/s - Next max velocity
        'linear',V_lim(end)),v_eng_max);
    F_grav(i) = m*g*sin(theta(i));                  %   N - Next gravitational force
    F_road(i) = m*g*Cr*cos(theta(i));               %   N - Next rolling resistance force
    F_aero(i) = k_D*V(i)^2;                         %   N - Next aerodynamic force
end

%%  Hard Constraints
if t_mode == 1
    c_ineq(1,:) = s_n-ds;                           %   Hard constraint on minimum following distance
    c_ineq(2,:) = (0-V).^3;
else
    c_ineq(1,:) = V-v_lim;                               %   Hard constraint on velocity
    c_ineq(2,:) = (0-V).^3;
end

c_eq = [];

end


function return_val = interpone(X_data, Y_data, X_q, str, extrapval)
    if X_q >= X_data(1) && X_q <= X_data(end)
        [~, I] = min((X_q - X_data).^2);
        if X_data(I) > X_q
            x_ub = X_data(I);
            x_lb = X_data(I-1);
            ind_xub = I;
            ind_xlb = I-1;
        else
            x_lb = X_data(I);
            ind_xlb = I;
            if I+1 >= length(X_data)    % M is probably equal to X_q and also the max value 
                x_lb = X_data(I-1);
                x_ub = I;
                ind_xub = I;
                ind_xlb = I-1;
            else
                x_ub = X_data(I+1);
                ind_xub = I+1;
            end
        end 
        y_0 = Y_data(ind_xlb);
        y_1 = Y_data(ind_xub);
        return_val = y_0 + (X_q - x_lb)*((y_1 - y_0)/(x_ub - x_lb));  
    else
        return_val = extrapval;
    end
end