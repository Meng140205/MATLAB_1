function plot_closed_loop_30cm(log,outputDir)
%PLOT_CLOSED_LOOP_30CM Save trajectory, tracking, and control plots.
if ~exist(outputDir,'dir'),mkdir(outputDir);end
t=log.time_s;
axesNames={'x','y','z'};
fig=figure('Name','Rod position tracking','Color','w','Visible','on');
for i=1:3
    subplot(3,1,i);
    if i<3
        scale=1000;offset=log.reference_m(1,i);unit='mm';
    else
        scale=100;offset=log.reference_m(1,i);unit='cm';
    end
    plot(t,scale*(log.reference_m(:,i)-offset),'k--','LineWidth',1.4);hold on
    plot(t,scale*(log.actual_m(:,i)-offset),'b','LineWidth',1.2);
    grid on;ylabel(sprintf('%s displacement (%s)',axesNames{i},unit));
    if i==1,legend('Reference','Simulated plant','Location','best');end
    ylim([-1,1]);
end
xlabel('Time (s)');
exportgraphics(fig,fullfile(outputDir,'01_rod_xyz_trajectory.png'),'Resolution',160);
savefig(fig,fullfile(outputDir,'01_rod_xyz_trajectory.fig'));%close(fig);

fig=figure('Name','Tracking errors','Color','w','Visible','on');
subplot(3,1,1);
plot(t,1000*log.error_m,'LineWidth',1.1);grid on
ylabel('Position error (mm)');legend('e_x','e_y','e_z','Location','best');
subplot(3,1,2);
plot(t,1000*vecnorm(log.error_m,2,2),'k','LineWidth',1.2);hold on
grid on
ylim([0 max(.2,1.2*max(1000*vecnorm(log.error_m,2,2)))]);
title('Zoomed error; Gazebo position guard is 3 mm');
ylabel('3D error (mm)');
subplot(3,1,3);
plot(t,log.angle_error_deg,'b','LineWidth',1.1);hold on
grid on
ylim([0 max(.05,1.2*max(log.angle_error_deg))]);
title('Zoomed attitude; Gazebo guard is 1 deg');
ylabel('Rod attitude (deg)');xlabel('Time (s)');
exportgraphics(fig,fullfile(outputDir,'02_tracking_errors.png'),'Resolution',160);
savefig(fig,fullfile(outputDir,'02_tracking_errors.fig'));%close(fig);

fig=figure('Name','Grasp and joint tracking','Color','w','Visible','on');
subplot(2,1,1);
plot(t,log.joint_error_rad,'b','LineWidth',1.1);hold on
grid on
ylim([0 max(.001,1.2*max(log.joint_error_rad))]);
title('Zoomed joint error; Gazebo guard is 0.02 rad');
ylabel('Max joint error (rad)');
subplot(2,1,2);
plot(t,log.grasp_error_mm,'b','LineWidth',1.1);grid on
ylabel('Max grasp position error (mm)');xlabel('Time (s)');
exportgraphics(fig,fullfile(outputDir,'03_grasp_and_joint.png'),'Resolution',160);
savefig(fig,fullfile(outputDir,'03_grasp_and_joint.fig'));%close(fig);

tc=t(1:numel(log.solve_ms));
fig=figure('Name','MPC effort and solve time','Color','w','Visible','on');
subplot(3,1,1);
plot(tc,100*log.torque_fraction,'b');hold on
yline(100,'r--','Torque bound');grid on
ylabel('Peak torque use (%)');
subplot(3,1,2);
plot(tc,log.torque_step_Nm,'b');hold on
yline(.9,'r--','0.9 N m / step');grid on
ylabel('Peak step (N m)');
subplot(3,1,3);
semilogy(tc,max(log.solve_ms,.1),'b');hold on
yline(50,'r--','50 ms period');grid on
ylabel('QP wall time (ms, log)');xlabel('Simulation time (s)');
title(sprintf('%d / %d solves exceeded 50 ms; max %.0f ms',...
    sum(log.solve_ms>50),numel(log.solve_ms),max(log.solve_ms)));
exportgraphics(fig,fullfile(outputDir,'04_control_and_solver.png'),'Resolution',160);
savefig(fig,fullfile(outputDir,'04_control_and_solver.fig'));%close(fig);
end
