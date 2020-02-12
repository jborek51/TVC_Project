%%  Plot Sim Comparison 
load('TVC_02-12_Results_C0_1.mat')
Sim0 = Sim;
load('TVC_02-12_Results_C1_1.mat')
Sim1 = Sim;
tF = 3.4;
figure()
subplot(2,2,1)
hold on;    grid on;
idx = find(Sim0.time == tF,1,'first');
plot(Sim0.Y_veh(1:idx),Sim0.Z_veh(1:idx),'r-')
idx = find(Sim1.time == tF,1,'first');
plot(Sim1.Y_veh(1:idx),Sim1.Z_veh(1:idx),'b-')
xlim([0 max([max(Sim.Y_veh) max(Sim0.Y_veh) max(Sim1.Y_veh)])]);
xlabel('Position [m]'); ylabel('Altitude [m]');
legend('None','SFB')
ax = gca;    ax.FontSize = 16;
subplot(2,2,3)
hold on;    grid on;
plot(Sim0.time,zeros(length(Sim0.time),1),'r-')
plot(Sim1.time,Sim1.SFB*180/pi,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('$\delta$ [$^o$]');
ax = gca;    ax.FontSize = 16;
subplot(2,2,2)
hold on;    grid on;
plot(Sim0.time,Sim0.phi_veh,'r-')
plot(Sim1.time,Sim1.phi_veh,'b-')
xlim([0 tF]);   ylim([-30 5])
xlabel('Time [s]'); ylabel('$\phi$ [$^o$]');
ax = gca;    ax.FontSize = 16;
subplot(2,2,4)
hold on;    grid on;
plot(Sim0.time,Sim0.p_veh,'r-')
plot(Sim1.time,Sim1.p_veh,'b-')
xlim([0 tF]);   ylim([-20 20])
xlabel('Time [s]'); ylabel('$p$ [$^o/\mathrm{s}$]');
ax = gca;    ax.FontSize = 16;