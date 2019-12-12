function MPC_S_Fun(block)

setup(block);

function setup(block)
%% Register number of input and output ports
block.NumInputPorts  = 4;
block.NumOutputPorts = 1;

%% Setup functional port properties to dynamically inherited.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Complexity   = 'Real';
block.InputPort(1).DatatypeID   = 0;
block.InputPort(1).Dimensions   = 7;
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Complexity   = 'Real';
block.InputPort(2).DatatypeID   = 0;  
block.InputPort(2).Dimensions   = 10;
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Complexity   = 'Real';
block.InputPort(3).DatatypeID   = 0;  
block.InputPort(3).Dimensions   = [28 2];
block.InputPort(3).DirectFeedthrough = true;

block.InputPort(4).Complexity   = 'Real';
block.InputPort(4).DatatypeID   = 0;  
block.InputPort(4).Dimensions   = 1;
block.InputPort(4).DirectFeedthrough = true;

block.OutputPort(1).Complexity   = 'Real';
block.OutputPort(1).SamplingMode   = 'Sample';
block.OutputPort(1).DataTypeId   = 0;
block.OutputPort(1).SamplingMode = 'Sample';
block.OutputPort(1).Dimensions   = 10;

%% Set block sample time to inherited
% block.SampleTimes = [-1 0];

%% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

%% Run accelerator on TLC
% block.SetAccelRunOnTLC(false);

%% Register methods
block.RegBlockMethod('SetInputPortSamplingMode', @sample);
block.RegBlockMethod('Outputs',                 @Output);  

%endfunction

function sample(block, idx, fd)
block.InputPort(idx).SamplingMode = fd;
block.OutputPort(1).SamplingMode  = fd;

function Output(block)
%%  Parameters
th_0 = block.InputPort(1).Data(1);                  %   rad - Current angle
thd_0 = block.InputPort(1).Data(2);                 %   rad/s - Current angular rate
v_0 = block.InputPort(1).Data(3);                     %   m/s - Body-y velocity 
w_0 = block.InputPort(1).Data(4);                     %   m/s - Body-z velocity 
% y_0 = block.InputPort(1).Data(5);                     %   m/s - Body-y velocity 
% z_0 = block.InputPort(1).Data(6);                     %   m/s - Body-z velocity 
t_0 = block.InputPort(1).Data(7);                   %   s - Sim time
param = block.InputPort(2).Data;                    %   Parameters 
F_th_v = block.InputPort(3).Data(:,1);                %   N - Thrust vector 
t_th = block.InputPort(3).Data(:,2);                %   s - Thrust time vector 
v_w = block.InputPort(4).Data;                      %   m/s - Wind velocity 

% g = param(1);                                       %   m/s^2 - Accel due to gravity
% rho = param(7);                                     %   kg/m^3 - Air density
% Cs = param(4);                                      %   Drag Coefficient
% A_side = param(5);                                  %   m^2 - Reference area
% Cd = param(2);                                      %   Drag Coefficient
% A_ref = param(3);                                  %   m^2 - Reference area
% k_S = .5*Cs*rho*A_side;                             %   kg/m - concentrated aero coef 
% k_D = .5*Cd*rho*A_ref;                             %   kg/m - concentrated aero coef 
% L1 = param(8);                                      %   m - Thrust moment arm 
% L2 = param(9);                                      %   m - Aero moment arm 
% J = param(6);                                       %   kg.m^2 - Polar moment of inertia
% m = param(10);                                       %   kg - Vehicle mass 

%%  Controller Settings
t_pred = 2;                                         %   s - Prediction horizon length
t_ctrl = 1;                                         %   s - Control horizon length
dt = .1;                                            %   s - Time step
np_steps = t_pred/dt;                               %   # of pred steps
nc_steps = t_ctrl/dt;                               %   # of control steps
% t_vec = t_0:dt:t_0+t_pred-dt;                       %   s - Horizon time vector 

%%  Initialize vectors
% th = zeros(1,np_steps);                             %   rad - Initialize theta vector
% thd = zeros(1,np_steps);                            %   rad/s - Initialize theta-dot vector
% th1 = zeros(1,np_steps);                            %   rad - Initialize desired theta vector
% thd1 = zeros(1,np_steps);                           %   rad/s - Initialize desired theta-dot vector
% v = zeros(1,np_steps);                              %   m/s - Initialize y-velocity vector
% w = zeros(1,np_steps);                              %   m/s - Initialize z-velocity vector
% y = zeros(1,np_steps);                              %   m/s - Initialize velocity vector
% z = zeros(1,np_steps);                              %   m/s - Initialize velocity vector
U = zeros(1,np_steps);                              %   rad - Initialize control vector
% F_th = zeros(1,np_steps);                           %   N - Initialize thrust force vector
% F_thy = zeros(1,np_steps);                          %   N - Initialize y-thrust force vector
% F_thz = zeros(1,np_steps);                          %   N - Initialize z-thrust force vector
% F_s = zeros(1,np_steps);                            %   N - Initialize side force vector       
% F_d = zeros(1,np_steps);                            %   N - Initialize drag force vector       
% F_gy = zeros(1,np_steps);                           %   N - Initialize y-gravity force vector       
% F_gz = zeros(1,np_steps);                           %   N - Initialize z-gravity force vector       

%%  Initial Parameters
% th(1) = th_0;   thd(1) = thd_0;                                     %   m - Initial position
% v(1) = v_0;     w(1) = w_0;
% % y(1) = y_0;     z(1) = x_0;
% F_th(1) = interp1(t_th,F_th_v,t_0,'linear','extrap'); %   N - Initial thrust force
% R = [cos(th(1)) -sin(th(1));sin(th(1)) cos(th(1))]; %   Ground to body rotation matrix 
% V = [v_0;w_0];  V_w = R'*[v_w;0];   V_app = V_w-V;
% if V_app(2) == 0
%     alpha = pi/2;
% elseif V_app(2) < 0
%     alpha = atan(V_app(1)/-V_app(2));
% else
%     alpha = pi+atan(V_app(1)/-V_app(2));
% end
% v_app = norm(V_app);
% v_ay = v_app*sin(alpha);
% v_az = v_app*cos(alpha);
% F_s(1) = k_S*v_ay^2*sign(v_ay);                             %   N - Initial aerodynamic force
% F_d(1) = k_D*v_az^2*sign(v_az);                             %   N - Initial aerodynamic force
% F_gy(1) = m*g*sin(th(1));
% F_gz(1) = m*g*cos(th(1));
% thd1(1) = 2/dt*(th1(1)-th(1))-thd(1); 
% U(1) = asin(((thd1(1)-thd(1))*J/dt...
%     -F_s(1)*L2)/(F_th(1)*L1));
% F_thy(1) = F_th(1)*sin(U(1));
% F_thz(1) = F_th(1)*cos(U(1));

%%  Loop through horizon for initial guess
% for i = 2:nc_steps                           
%     v(i) = v(i-1)+dt/m*(-m*w(i-1)*thd(i-1)...
%         +(F_thy(i-1)+F_s(i-1)+F_gy(i-1)));
%     w(i) = w(i-1)+dt/m*(m*v(i-1)*thd(i-1)...
%         +(F_thz(i-1)-F_d(i-1)-F_gz(i-1)));
%     thd(i) = thd(i-1)+dt/J*(F_thy(i-1)*L1...
%         +F_s(i-1)*L2);
%     th(i) = th(i-1)+(thd(i)+thd(i-1))/2*dt;
%     F_th(i) = interp1(t_th,F_th_v,t_vec(i),...
%         'linear','extrap');
%     R = [cos(th(1)) -sin(th(1));sin(th(1)) cos(th(1))];
%     V = [v(i);w(i)];  V_w = R'*[v_w;0];   V_app = V_w-V;
%     if V_app(2) == 0
%         alpha = pi/2;
%     elseif V_app(2) < 0
%         alpha = atan(V_app(1)/-V_app(2));
%     else
%         alpha = pi+atan(V_app(1)/-V_app(2));
%     end
%     v_app = norm(V_app);
%     v_ay = v_app*sin(alpha);
%     v_az = v_app*cos(alpha);
%     F_s(i) = k_S*v_ay^2*sign(v_ay);                             %   N - Initial aerodynamic force
%     F_d(i) = k_D*v_az^2*sign(v_az);                             %   N - Initial aerodynamic force
%     F_gy(i) = m*g*sin(th(1));
%     F_gz(i) = m*g*cos(th(1));
%     thd1(i) = 2/dt*(th1(i)-th(i))-thd(i); 
%     U(i) = asin(((thd1(i)-thd(i))*J/dt...
%         -F_s(i)*L2)/(F_th(i)*L1));
%     F_thy(i) = F_th(i)*sin(U(i));
%     F_thz(i) = F_th(i)*cos(U(i));
% end

%%  Set up and perform Optimization                 
states = [th_0,thd_0,v_0,w_0,t_0];
thrust = [F_th_v,t_th];
lb = -2*pi/180.*ones(1,nc_steps);                  %   N - Control lower bound
ub = 2*pi/180.*ones(1,nc_steps);                   %   N - Control upper bound
J = @(U_0)obj_fun(U_0,states,param,thrust,v_w);     %   Initialize objective function 
options = optimoptions(@fmincon,'MaxIterations'...  %   Fmincon options
    ,500,'MaxFunctionEvaluations',5000000,...
    'Algorithm','sqp','Display','iter');

U_new = fmincon(J,U(1:nc_steps),...                 %   rad - Optimal control sequence 
    [],[],[],[],lb,ub,[],options);

% U_app = zeros(1,np_steps);                          %   rad - Initialize optimal control sequence 
% U_app(1:nc_steps) = U_new;                          %   rad - Optimal control sequence 
% U_app(nc_steps+1:np_steps) = U_app(nc_steps)...     %   rad - Optimal control sequence 
%     .*ones(1,np_steps-nc_steps);

%%  Output Variables
block.OutputPort(1).Data = U_new;                   %   rad - New control sequence 

%endfunction

