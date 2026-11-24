function report = run_closed_loop_30cm(varargin)
%RUN_CLOSED_LOOP_30CM Offline 30 cm dual-arm/rod closed-loop MPC simulation.
% Uses the stage 4.3/4.4 mathematical model, not ROS or Gazebo.
% The plant state is integrated from commanded torques; it is never set to
% the reference trajectory. See README.md for the plant-model limitation.

p = inputParser;
addParameter(p,'outputDir',fullfile(fileparts(mfilename('fullpath')),'results'));
addParameter(p,'simDuration',62,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'disturbanceNm',0.02,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'makePlots',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'verbose',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt=p.Results;
root=fileparts(mfilename('fullpath'));
if ~exist(opt.outputDir,'dir'),mkdir(opt.outputDir);end

urdf=fullfile(root,'model_data','dualarm_effort.urdf');
sdf=fullfile(root,'model_data','rod_with_pose.sdf');
arms={dualarm.loadArm(urdf,'left',[.01 .01]),...
      dualarm.loadArm(urdf,'right',[.01 .01])};
payload=dualarm.rodInertia(sdf);
q0=[1.1519855608624812;-2.1271652434076493;...
    2.0510471815526485;-1.4927644876148773;...
    1.570789063937433;-0.41880370305441245;...
    1.1519852783099607;-2.127165249317388;...
    2.0510485070664495;-1.492764499929628;...
    1.5707890042603168;-0.41880370006208273];
rod0=eye(4);
rod0(1:3,1:3)=dualarm.quatToRot(...
    [-0.0008759998879626704;9.99999616311858e-7;...
     8.759998879629624e-10;0.9999996163115246]);
rod0(1:3,4)=[2e-6;0.299657;1.27105];

height=.30;duration=60;dt=.05;N=10;
path=dualarm.liftReference(arms,q0,rod0,height,duration);
grasp=cell(1,2);
for i=1:2
    sel=(6*i-5):(6*i);
    [tip,~]=dualarm.armGeometry(arms{i},q0(sel));
    grasp{i}=rod0\tip;
end
lower=[arms{1}.lower;arms{2}.lower];
upper=[arms{1}.upper;arms{2}.upper];
limits=[arms{1}.effort;arms{2}.effort];
mu=[arms{1}.friction;arms{2}.friction];

q=q0;dq=zeros(12,1);rod=rod0;rodVelocity=zeros(6,1);
initial=dualarm.localModel(arms,q,dq,rod,rodVelocity,...
    grasp,payload,zeros(12,1));
previous=-pinv(initial.F)*initial.g;
equilibriumResidual=norm(initial.F*previous+initial.g);
if equilibriumResidual>1e-3 || any(abs(previous)>limits)
    error('Initial grasp is not a feasible static equilibrium');
end
biasState=[];

steps=round(opt.simDuration/dt);
time=(0:steps)'*dt;
actual=zeros(steps+1,3);desired=actual;errorXYZ=actual;
qActual=zeros(steps+1,12);qDesired=qActual;
angleDeg=zeros(steps+1,1);graspMm=angleDeg;
jointError=angleDeg;rodSpeed=angleDeg;
torque=zeros(steps,12);torqueFraction=zeros(steps,1);
torqueStep=torqueFraction;solveMs=torqueFraction;
qpResidual=torqueFraction;biasNorm=torqueFraction;
method=cell(steps,1);last=0;stopReason='';
for k=1:steps
    t=time(k);
    [qref,dqref,ref,refVelocity]=dualarm.evalLiftReference(path,t);
    actual(k,:)=rod(1:3,4).';
    desired(k,:)=ref(1:3,4).';
    errorXYZ(k,:)=actual(k,:)-desired(k,:);
    qActual(k,:)=q.';qDesired(k,:)=qref.';
    angleDeg(k)=rad2deg(norm(dualarm.rotLog(...
        rod(1:3,1:3)*ref(1:3,1:3).')));
    jointError(k)=max(abs(q-qref));
    rodSpeed(k)=norm(rodVelocity(1:3));
    [biasState,bias]=dualarm.trackingBiasStep(...
        biasState,qref-q,t,t>0 && t<duration);
    [friction,~]=dualarm.frictionFeedforward(mu,dqref,dq);
    try
        prediction=dualarm.localModel(arms,q,dq,rod,rodVelocity,...
            grasp,payload,friction+bias);
        refs=zeros(N+1,36);
        for h=0:N
            [qr,dqr,pr,vr]=dualarm.evalLiftReference(path,t+h*dt);
            poseError=[pr(1:3,4)-rod(1:3,4);...
                dualarm.rotLog(pr(1:3,1:3)*rod(1:3,1:3).')];
            refs(h+1,:)=[qr-q;poseError;dqr;vr].';
        end
        solved=dualarm.solveMPC(prediction,q,refs,previous,...
            lower,upper,limits,N,dt);
        u=solved.torque;
        % Separate plant evaluation: a small, fixed, bounded unmodelled
        % joint disturbance prevents an exact model-equals-plant replay.
        disturbance=opt.disturbanceNm*sin(.13*t+(1:12)'*.4);
        plant=dualarm.localModel(arms,q,dq,rod,rodVelocity,...
            grasp,payload,friction+disturbance);
        acceleration=plant.F*u+plant.g;
        [q,dq,rod,rodVelocity]=advancePlant(...
            q,dq,rod,rodVelocity,acceleration,dt);
        if any(~isfinite([q;dq;rod(:);rodVelocity]))
            error('Plant state became non-finite');
        end
        graspMm(k)=graspErrorMm(arms,q,rod,grasp);
        torque(k,:)=u.';
        torqueFraction(k)=max(abs(u)./limits);
        torqueStep(k)=max(abs(u-previous));
        solveMs(k)=1000*solved.elapsed;
        qpResidual(k)=solved.residual;
        biasNorm(k)=max(abs(bias));
        method{k}=solved.method;
        previous=u;last=k;
    catch ME
        stopReason=sprintf('t=%.2f s: %s',t,ME.message);
        warning('Closed-loop simulation stopped: %s',stopReason);
        last=k-1;
        break
    end
    if opt.verbose && (mod(k,100)==0 || k==steps)
        fprintf('t=%5.1f s  lift=%7.2f mm  |e_xyz|=%6.3f mm  max|dq|=%6.4f rad/s\n',...
            t,1000*(rod(3,4)-rod0(3,4)),...
            1000*norm(errorXYZ(k,:)),max(abs(dq)));
    end
end

finalIndex=last+1;
if finalIndex<=numel(time)
    t=time(finalIndex);
    [qref,~,ref,~]=dualarm.evalLiftReference(path,t);
    actual(finalIndex,:)=rod(1:3,4).';
    desired(finalIndex,:)=ref(1:3,4).';
    errorXYZ(finalIndex,:)=actual(finalIndex,:)-desired(finalIndex,:);
    qActual(finalIndex,:)=q.';qDesired(finalIndex,:)=qref.';
    angleDeg(finalIndex)=rad2deg(norm(dualarm.rotLog(...
        rod(1:3,1:3)*ref(1:3,1:3).')));
    jointError(finalIndex)=max(abs(q-qref));
    rodSpeed(finalIndex)=norm(rodVelocity(1:3));
end

valid=1:finalIndex;
log=struct('time_s',time(valid),'actual_m',actual(valid,:),...
    'reference_m',desired(valid,:),'error_m',errorXYZ(valid,:),...
    'q_rad',qActual(valid,:),'q_ref_rad',qDesired(valid,:),...
    'angle_error_deg',angleDeg(valid),'grasp_error_mm',graspMm(valid),...
    'joint_error_rad',jointError(valid),'rod_speed_m_s',rodSpeed(valid),...
    'torque_Nm',torque(1:last,:),'torque_fraction',torqueFraction(1:last),...
    'torque_step_Nm',torqueStep(1:last),'solve_ms',solveMs(1:last),...
    'qp_residual',qpResidual(1:last),'bias_Nm',biasNorm(1:last),...
    'qp_method',{method(1:last)});
matPath=fullfile(opt.outputDir,'closed_loop_30cm.mat');
save(matPath,'log','path','grasp','q0','rod0','equilibriumResidual');
csvPath=fullfile(opt.outputDir,'closed_loop_30cm.csv');
tab=table(log.time_s,log.actual_m(:,1),log.actual_m(:,2),...
    log.actual_m(:,3),log.reference_m(:,1),log.reference_m(:,2),...
    log.reference_m(:,3),log.error_m(:,1),log.error_m(:,2),...
    log.error_m(:,3),log.joint_error_rad,log.angle_error_deg,...
    log.rod_speed_m_s,'VariableNames',{'time_s','x_m','y_m','z_m',...
    'x_ref_m','y_ref_m','z_ref_m','ex_m','ey_m','ez_m',...
    'max_joint_error_rad','angle_error_deg','rod_speed_m_s'});
writetable(tab,csvPath);
if opt.makePlots
    plot_closed_loop_30cm(log,opt.outputDir);
end
maxError=1000*max(vecnorm(log.error_m,2,2));
finalLift=rod(3,4)-rod0(3,4);
completed=last==steps;
accepted=completed && abs(finalLift-height)<=.002 && ...
    maxError<=3 && max(log.joint_error_rad)<=.02 && ...
    max(log.angle_error_deg)<=1 && max(abs(dq))<=.01 && ...
    max(abs(rodVelocity))<=.005;
report=struct('completed',completed,'offline_acceptance',accepted,...
    'simulated_until_s',time(finalIndex),'desired_lift_m',height,...
    'actual_lift_m',finalLift,'max_xyz_error_mm',maxError,...
    'max_joint_error_rad',max(log.joint_error_rad),...
    'max_angle_error_deg',max(log.angle_error_deg),...
    'max_torque_fraction',max([0;log.torque_fraction]),...
    'max_torque_step_Nm',max([0;log.torque_step_Nm]),...
    'max_solver_ms',max([0;log.solve_ms]),...
    'median_solver_ms',median(log.solve_ms),...
    'p95_solver_ms',prctile(log.solve_ms,95),...
    'solver_over_50ms',sum(log.solve_ms>50),...
    'max_grasp_error_mm',max(log.grasp_error_mm),...
    'equilibrium_residual_m_s2',equilibriumResidual,...
    'stop_reason',stopReason,'mat_file',matPath,'csv_file',csvPath);
disp(report);
end

function [q,dq,rod,twist]=advancePlant(q,dq,rod,twist,a,h)
q=q+h*dq+.5*h*h*a(1:12);
rod(1:3,4)=rod(1:3,4)+h*twist(1:3)+.5*h*h*a(13:15);
angle=h*twist(4:6)+.5*h*h*a(16:18);
theta=norm(angle);
K=dualarm.skew(angle);
if theta<1e-8
    R=eye(3)+K+.5*K*K;
else
    R=eye(3)+(sin(theta)/theta)*K+...
        ((1-cos(theta))/(theta*theta))*K*K;
end
rod(1:3,1:3)=R*rod(1:3,1:3);
dq=dq+h*a(1:12);
twist=twist+h*a(13:18);
end

function value=graspErrorMm(arms,q,rod,grasp)
value=0;
for i=1:2
    sel=(6*i-5):(6*i);
    [tip,~]=dualarm.armGeometry(arms{i},q(sel));
    expected=rod*grasp{i};
    value=max(value,1000*norm(tip(1:3,4)-expected(1:3,4)));
end
end
