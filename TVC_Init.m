%%  Thrust-vector Control Project Initialization File
%   Created by John Borek
%   
%   This flight simulation has 3 degrees of freedom: vertical (z),
%   horizontal (y), and rotational (theta). This file initializes all
%   parameters, runs a simulation, saves the outputs, and plots the
%   results.
clear;  clc

%%   Environmental Parameters
ENV.g = 9.80665;                                %   m/s^2 - Gravitational acceleration
ENV.rho_b = 1.225;                              %   kg/m^3 - Air density at sea level
ENV.mu = 1.983e-5;                              %   Pa.s - Dynamic velocity of air
ENV.nu_a = 2.503e-5;                            %   Coefficient to kinetic viscosity
ENV.P_b = 101325;                               %   Pa - Air Pressure at sea level
ENV.T = 288.15;                                 %   K - Air temperature
ENV.M = 0.0289644;                              %   kg/mol - Molar mass of air
ENV.R = 8.31432;                                %   N.m/(mol.K) - Universal gas constant for air
ENV.v_w0 = 20*.447;                             %   mph - Wind velocity at ground level
ENV.v_wf = 20*.447;                             %   mph - Wind velocity at 5000 ft
ENV.z_0 = 0;                                    %   m - Initial altitude above sea level
ENV.z_r = 36*.0254;                             %   in -> m - Rail length
ENV.theta_0 = 0*pi/180;                         %   deg -> rad - Launch angle off of vertical
ENV.t_f = 50;                                   %   s - Simulation end time

%%  Load Motor Parameters
load('F15.mat');                                %   Estes F15 motor
Motor.k_m = 1.0;                                %   Motor variance constant

%%   Rocket Parameters
VEH.m_lb = 2;                                   %   lb - Launch vehicle's unloaded mass
VEH.m = VEH.m_lb*0.453592;                      %   kg - Launch vehicle's unloaded mass
VEH.C_d = .46;                                  %   Coefficient of drag
VEH.C_s = .5;                                   %   Coefficient of side force
VEH.L = 37*.0254;                               %   in -> m - Rocket length
VEH.D = 3.1*.0254;                              %   in -> m - Major airframe diameter
VEH.A_ref = pi*VEH.D^2/4;                       %   m^2 - Drag reference area
VEH.A_side = VEH.L*VEH.D;                       %   m^2 - Side force reference area
VEH.CG_1 = 21.70*.0254;                         %   in -> m - Center of gravity on pad
VEH.CG_2 = 17.60*.0254;                         %   in -> m - Center of gravity after burnout
VEH.L2 = VEH.CG_1 - VEH.CG_2;                   %   m - Change in CG position 
VEH.CP = 17.78*.0254;                           %   in -> m - Center of Pressure
VEH.J = .0625*VEH.m*VEH.D^2+1/12*VEH.m*VEH.L^2; %   kg.m^2 - Polar moment of inertia

%%  Actuator Parameters
ACT.u_max = 20*pi/180;                          %   deg -> rad - Max servo angle
ACT.u_min = -20*pi/180;                         %   deg -> rad - Min servo angle
ACT.w_s = 30*pi/180;                            %   deg/s -> rad/s - Servo rate limit
ACT.tau = .1;                                   %   s - Actuator time delay

%%   Sensor Parameters
SENS.tau_alt = .004;                            %   s - Altimeter time delay
SENS.tau_accel = .004;                          %   s - Accelerometer time delay
SENS.kn_alt = .01;                              %   Altimeter noise constant
SENS.kn_accel = .05;                            %   Accelerometer noise constant

%%  Controller Settings
CTRL.k_s = 2;                                   %   Control switch:  1 = off;  2 = SFB;  3 = MPC
CTRL.Kp = .25;                                  %   rad/(m/s) - Proportional gain
CTRL.Ki = .1;                                   %   rad/m - Integral gain
CTRL.Kd = 0;                                    %   rad/(m/s^2) - Derivative gain
CTRL.tau = .1;                                  %   Filter time constant
CTRL.K1 = 1;                                    %   rad/rad - Feedback gain on orientation 
CTRL.K2 = 1;                                    %   rad/(rad/s) - Feedback gain on rotation rate

%%  Linearized Model
LIN.A = [-6.94*VEH.C_s*VEH.A_side*ENV.rho_b/VEH.m 0 -8 0 0 0
    0 -8*VEH.C_d*VEH.A_ref*ENV.rho_b/VEH.m 2 0 0 0
    -6.94*VEH.C_s*VEH.A_side*ENV.rho_b*(-VEH.CP+mean([VEH.CG_1 VEH.CG_2]))/VEH.J 0 0 0 0 0
    1 0 0 0 0 -8;0 1 0 0 0 2;0 0 1 0 0 0];
LIN.B = [13.84/VEH.m 0 6.94*VEH.C_s*VEH.A_side*ENV.rho_b/VEH.m
    0 1/VEH.m 0
    13.84*(VEH.L-mean([VEH.CG_1 VEH.CG_2]))/VEH.J 0 6.94*VEH.C_s*VEH.A_side*ENV.rho_b*(-VEH.CP+mean([VEH.CG_1 VEH.CG_2]))/VEH.J
    0 0 0;0 0 0;0 0 0];
LIN.C = [0 0 1 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
LIN.D = zeros(3);

%%  Run Simulation
sim('TVC_Model');

%%  Save Results
Computer = 1;
run Log_Results
run getFilename
save(strcat(parent,'\',filename),'Sim')

%%  Plots
run Plot_Results