%%
% tFinal = Sim.time(end);
tF = 4;
V_W = round((ENV.v_w0+ENV.v_wf)/2/.447);
% idx = find(Sim.time == tF,1,'first');
idx = length(Sim.time);
%%  Plot translational states
figure()
subplot(3,3,2)
hold on;    grid on;
% plot(Sim.time,Sim.Z_veh,'c-')
% plot(Sim.time,Sim.Y_veh,'m-')
% xlim([0 tF]);
% xlabel('Time [s]'); ylabel('Distance [m]');
plot(Sim.Y_veh(1:idx),Sim.Z_veh(1:idx),'k-')
xlim([0 max(Sim.Y_veh)]);
xlabel('Position [m]'); ylabel('Altitude [m]');
ax = gca;    ax.FontSize = 12;
% text(20,10,'$v_w\longrightarrow$','FontSize',16)
subplot(3,3,5)
hold on;    grid on;
plot(Sim.time,Sim.Y_dot_b,'r-')
plot(Sim.time,Sim.Z_dot_b,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
% legend('Z','Y');
ax = gca;    ax.FontSize = 12;
subplot(3,3,8)
hold on;    grid on;
plot(Sim.time,Sim.Y_ddot_b,'r-')
plot(Sim.time,Sim.Z_ddot_b,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Body Accel [m/$\mathrm{s^2}$]');
% L = legend('$Z$','$Y$');
ax = gca;    ax.FontSize = 12;

%%  Plot rotational states 
% figure()
subplot(3,3,3)
hold on;    grid on;
plot(Sim.time,Sim.theta_veh,'k-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\theta$ [$^o$]');
ax = gca;    ax.FontSize = 12;
subplot(3,3,6)
hold on;    grid on;
plot(Sim.time,Sim.theta_dot_veh,'k-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\omega$ [$^o/\mathrm{s}$]');
ax = gca;    ax.FontSize = 12;
subplot(3,3,9)
hold on;    grid on;
plot(Sim.time,Sim.theta_ddot_veh,'k-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\dot{\omega}$ [$^o/\mathrm{s}^2$]');
ax = gca;    ax.FontSize = 12;
%%  Plot Forces
% figure()
subplot(3,3,1)
hold on;    grid on;
plot(Sim.time,Sim.Y_thrust,'r-')
plot(Sim.time,Sim.Z_thrust,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Thrust [N]');
legend('Y','Z');
ax = gca;    ax.FontSize = 12;
subplot(3,3,4)
hold on;    grid on;
plot(Sim.time,Sim.F_side,'r-')
plot(Sim.time,Sim.Y_grav,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Force [N]');
legend('$F_\mathrm{side}$','$F_\mathrm{grav}^y$','location','northwest');
ax = gca;    ax.FontSize = 12;
subplot(3,3,7)
hold on;    grid on;
plot(Sim.time,Sim.Y_thrust+Sim.F_side-Sim.Y_grav,'r-')
plot(Sim.time,Sim.Z_thrust-Sim.F_drag-Sim.Z_grav+Sim.Norm,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\sum F$ [N]');
ax = gca;    ax.FontSize = 12;
% T = suptitle(['TVC Simulation: $m$ = ',num2str(VEH.m_lb),'kg, $v_w$ = ',num2str(V_W),'mph, F15 Motor']);
% % T.FontSize = 14;    
% T.Position = [.5 -.075 0];
