function [J] = obj_fun(U_0,states,param,thrust,v_w)
%	MPC objective function for braking penalty

th_0 = states(1);                                   %   rad - Current angle
thd_0 = states(2);                                  %   rad/s - Current angular rate
v_0 = states(3);                                    %   m/s - Body-y velocity 
w_0 = states(4);                                    %   m/s - Body-z velocity 
t_0 = states(5);                                    %   s - Sim time
F_th_v = thrust(:,1);                               %   N - Thrust vector 
t_th = thrust(:,2);                                 %   s - Thrust time vector 

g = param(1);                                       %   m/s^2 - Accel due to gravity
rho = param(7);                                     %   kg/m^3 - Air density
Cs = param(4);                                      %   Side Coefficient
A_side = param(5);                                  %   m^2 - Reference area
Cd = param(2);                                      %   Drag Coefficient
A_ref = param(3);                                   %   m^2 - Reference area
k_S = .5*Cs*rho*A_side;                             %   kg/m - concentrated side coef 
k_D = .5*Cd*rho*A_ref;                              %   kg/m - concentrated drag coef 
L1 = param(8);                                      %   m - Thrust moment arm 
L2 = param(9);                                      %   m - Aero moment arm 
J = param(6);                                       %   kg.m^2 - Polar moment of inertia
m = param(10);                                      %   kg - Vehicle mass 

%%  Controller Settings
t_pred = 2;                                         %   s - Prediction horizon length
t_ctrl = 1;                                         %   s - Control horizon length
dt = .1;                                            %   s - Time step
np_steps = t_pred/dt;                               %   # of pred steps
nc_steps = t_ctrl/dt;                               %   # of control steps
t_vec = t_0:dt:t_0+t_pred-dt;                       %   s - Horizon time vector 

%%  Initialize vectors
th = zeros(1,np_steps);                             %   rad - Initialize theta vector
thd = zeros(1,np_steps);                            %   rad/s - Initialize theta-dot vector
v = zeros(1,np_steps);                              %   m/s - Initialize y-velocity vector
w = zeros(1,np_steps);                              %   m/s - Initialize z-velocity vector
U = zeros(1,np_steps);                              %   rad - Initialize control vector
U(1:nc_steps) = U_0;                                %   rad - Control vector
U(nc_steps+1:np_steps) = U_0(end)...                %   rad - Control vector 
    *ones(1,np_steps-nc_steps);
F_th = zeros(1,np_steps);                           %   N - Initialize thrust force vector
F_thy = zeros(1,np_steps);                          %   N - Initialize y-thrust force vector
F_thz = zeros(1,np_steps);                          %   N - Initialize z-thrust force vector
F_s = zeros(1,np_steps);                            %   N - Initialize side force vector       
F_d = zeros(1,np_steps);                            %   N - Initialize drag force vector       
F_gy = zeros(1,np_steps);                           %   N - Initialize y-gravity force vector       
F_gz = zeros(1,np_steps);                           %   N - Initialize z-gravity force vector       

%%  Initial Parameters
th(1) = th_0;   thd(1) = thd_0;                         %   Initial angular states
v(1) = v_0;     w(1) = w_0;                             %   Initial velocity states 
F_th(1) = interpone(t_th,F_th_v,t_0,'linear',0);        %   N - Initial thrust force
R = [cos(th(1)) -sin(th(1));sin(th(1)) cos(th(1))];     %   Ground to body rotation matrix 
V = [v(1);w(2)];  V_w = R'*[v_w;0];   V_app = V_w-V;    %   m/s - Wind vectors
if V_app(2) == 0
    alpha = pi/2;
elseif V_app(2) < 0
    alpha = atan(V_app(1)/-V_app(2));
else
    alpha = pi+atan(V_app(1)/-V_app(2));
end
v_app = norm(V_app);                                    %   m/s - Apparent wind magnitude 
v_ay = v_app*sin(alpha);                                %   m/s - Y-apparent wind
v_az = v_app*cos(alpha);                                %   m/s - Z-apparent wind
F_s(1) = k_S*v_ay^2*sign(v_ay);                         %   N - Initial side force
F_d(1) = k_D*v_az^2*sign(v_az);                         %   N - Initial drag force
F_gy(1) = m*g*sin(th(1));                               %   N - Initial y-grav force
F_gz(1) = m*g*cos(th(1));                               %   N - Initial z-grav force
F_thy(1) = F_th(1)*sin(U(1));                           %   N - Initial y-thrust force 
F_thz(1) = F_th(1)*cos(U(1));                           %   N - Initial z-thrust force 

%%  Loop through horizon for initial guess
for i = 2:nc_steps                           
    v(i) = v(i-1)+dt/m*(-m*w(i-1)*thd(i-1)...
        +(F_thy(i-1)+F_s(i-1)+F_gy(i-1)));
    w(i) = w(i-1)+dt/m*(m*v(i-1)*thd(i-1)...
        +(F_thz(i-1)-F_d(i-1)-F_gz(i-1)));
    thd(i) = thd(i-1)+dt/J*(F_thy(i-1)*L1...
        +F_s(i-1)*L2);
    th(i) = th(i-1)+(thd(i)+thd(i-1))/2*dt;
    F_th(i) = interpone(t_th,F_th_v,t_vec(i),...
        'linear',0);
    R = [cos(th(1)) -sin(th(1));sin(th(1)) cos(th(1))];
    V = [v(i);w(i)];  V_w = R'*[v_w;0];   V_app = V_w-V;
    if V_app(2) == 0
        alpha = pi/2;
    elseif V_app(2) < 0
        alpha = atan(V_app(1)/-V_app(2));
    else
        alpha = pi+atan(V_app(1)/-V_app(2));
    end
    v_app = norm(V_app);
    v_ay = v_app*sin(alpha);
    v_az = v_app*cos(alpha);
    F_s(i) = k_S*v_ay^2*sign(v_ay);                             %   N - Initial aerodynamic force
    F_d(i) = k_D*v_az^2*sign(v_az);                             %   N - Initial aerodynamic force
    F_gy(i) = m*g*sin(th(1));
    F_gz(i) = m*g*cos(th(1));
    F_thy(i) = F_th(i)*sin(U(i));
    F_thz(i) = F_th(i)*cos(U(i));
end
%%  Objective Function
J = 0*sum(U)^2+1e2*sum(th)^2+0*sum(thd)^2;

end

function return_val = interpone(X_data, Y_data, X_q, ~, extrapval)
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
                x_lb = X_data(I);
                x_ub = I+1;
                ind_xub = I+1;
                ind_xlb = I;
            else
                x_ub = X_data(I+1);
                ind_xub = I+1;
            end
        end 
        y_0 = Y_data(ind_xlb);
        y_1 = Y_data(ind_xub);
        return_val = y_0 + (X_q - x_lb)*((y_1 - y_0)/(x_ub - x_lb));  
    elseif ischar(extrapval) 
        return_val = Y_data(end)+(X_q-X_data(end))/(X_data(end-1)-X_data(end))*(Y_data(end-1)-Y_data(end));
    else
        return_val = extrapval;
    end
end
