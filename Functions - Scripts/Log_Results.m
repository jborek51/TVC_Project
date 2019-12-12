%%  Log results
Sim.time = VEH_out.time;
Sim.Z_veh = VEH_out.signals.values(:,6);
Sim.Z_dot_veh = VEH_out.signals.values(:,8);
% Sim.Z_ddot_veh = VEH_out.signals.values(:,3);
Sim.Y_veh = VEH_out.signals.values(:,5);
Sim.Y_dot_veh = VEH_out.signals.values(:,7);
% Sim.Y_ddot_veh = VEH_out.signals.values(:,6);
Sim.theta_veh = VEH_out.signals.values(:,11)*180/pi;
Sim.theta_dot_veh = VEH_out.signals.values(:,10)*180/pi;
Sim.theta_ddot_veh = VEH_out.signals.values(:,9)*180/pi;

% Sim.Z_b = VEH_out.signals.values(:,10);
Sim.Z_dot_b = VEH_out.signals.values(:,2);
Sim.Z_ddot_b = VEH_out.signals.values(:,4);
% Sim.Y_b = VEH_out.signals.values(:,13);
Sim.Y_dot_b = VEH_out.signals.values(:,1);
Sim.Y_ddot_b = VEH_out.signals.values(:,3);

Sim.F_thrust = Dist_out.signals.values(:,1);
Sim.Z_thrust = Dist_out.signals.values(:,2);
Sim.Y_thrust = Dist_out.signals.values(:,3);
Sim.T_thrust = Dist_out.signals.values(:,4);
Sim.m_motor = Dist_out.signals.values(:,5);
Sim.F_drag = Dist_out.signals.values(:,6);
Sim.F_side = Dist_out.signals.values(:,7);
Sim.T_aero = Dist_out.signals.values(:,8);
Sim.F_grav = Dist_out.signals.values(:,9);
Sim.Z_grav = Dist_out.signals.values(:,10);
Sim.Y_grav = Dist_out.signals.values(:,11);
Sim.Norm = Dist_out.signals.values(:,12);

Sim.delta = ACT_out.signals.values(:,1);
Sim.SFB = CTRL_out.signals.values(:,1);
Sim.MPC = CTRL_out.signals.values(:,2);
Sim.Sat_diff = ACT_out.signals.values(:,2);
Sim.V_app = V_app;
Sim.V_a = V_a;
Sim.V_appD = V_appD;
Sim.V_appS = V_appS;
Sim.alpha = alpha;
Sim.SM = Stability;
Sim.mass = VEH.m_lb;
Sim.date = datestr(now, 'dd-mmmm_HH-MM');
