%%
% tFinal = Sim.time(end);
tF = 4;
V_W = round((ENV.v_w0+ENV.v_wf)/2/.447);
idx = find(Sim.time == tF,1,'first');
%%  Plot translational states
figure()
subplot(3,3,2)
hold on;    grid on;
% plot(Sim.Y_veh(1:idx),Sim.Z_veh(1:idx),'k-')
plot(Sim.time,Sim.Z_veh,'c-')
plot(Sim.time,Sim.Y_veh,'m-')
% xlim([0 max(Sim.Y_veh)]);
xlim([0 tF]);
% xlabel('Position [m]'); ylabel('Altitude [m]');
xlabel('Time [s]'); ylabel('Distance [m]');
ax = gca;    ax.FontSize = 12;
% text(20,10,'$v_w\longrightarrow$','FontSize',16)
subplot(3,3,5)
hold on;    grid on;
plot(Sim.time,Sim.Z_dot_b,'c-')
plot(Sim.time,Sim.Y_dot_b,'m-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
% legend('Z','Y');
ax = gca;    ax.FontSize = 12;
subplot(3,3,8)
hold on;    grid on;
plot(Sim.time,Sim.Z_ddot_b,'c-')
plot(Sim.time,Sim.Y_ddot_b,'m-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Acceleration [m/$\mathrm{s^2}$]');
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
plot(Sim.time,Sim.F_thrust,'k-')
plot(Sim.time,Sim.Z_thrust,'c--')
plot(Sim.time,Sim.Y_thrust,'m--')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Thrust [N]');
legend('Nom','Z','Y');
ax = gca;    ax.FontSize = 12;
subplot(3,3,4)
hold on;    grid on;
plot(Sim.time,Sim.F_drag,'c-')
plot(Sim.time,Sim.F_side,'m-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Aero [N]');
ax = gca;    ax.FontSize = 12;
subplot(3,3,7)
hold on;    grid on;
plot(Sim.time,Sim.Z_thrust-Sim.F_drag-Sim.Z_grav+Sim.Norm,'c-')
% plot(Sim.time,Sim.Z_thrust-Sim.F_drag-Sim.Z_grav,'c-')
plot(Sim.time,Sim.Y_thrust+Sim.F_side-Sim.Y_grav,'m-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\sum F$ [N]');
ax = gca;    ax.FontSize = 12;
% T = suptitle(['TVC Simulation: $m$ = ',num2str(VEH.m_lb),'kg, $v_w$ = ',num2str(V_W),'mph, F15 Motor']);
% % T.FontSize = 14;    
% T.Position = [.5 -.075 0];
