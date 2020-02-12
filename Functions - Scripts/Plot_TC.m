%%  Plot Thrust Curve
figure()
hold on;    grid on;
plot(Sim.time,Sim.F_thrust,'b-')
xlim([0 tF]);
xlabel('Time [s]'); ylabel('Thrust [N]');
ax = gca;    ax.FontSize = 16;