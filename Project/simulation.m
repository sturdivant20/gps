%% GPS Project | Daniel Sturdivant
clc; clear; close all;

c = 299792458;

runKalman = true;
% runKalman = false;
% kalmanType = 'tight';
kalmanType = 'loose';

%% DATA LOADING

load("+data/simPerfect.mat");
psi_init = ref.psi_n2b(:,1);
user_lla = ecef2lla(ref.ecef');
user_vel = ref.ecef_vel;
user_xyz = ref.ecef;
t = ref.time;


%% INITIALIZATION

% imu bias
b_as = zeros(3,1);
b_gs = zeros(3,1);
b_ad = zeros(3,1);
b_gd = zeros(3,1);

% preallocation (for speed :])
Li = length(imu.time);
Lg = length(gps.gpsTime);
C = zeros(3,3,Li);
att = zeros(3,Li);
vel = zeros(3,Li);
pos = zeros(3,Li);
lla = zeros(3,Li);
att_mech = zeros(3,Li);
vel_mech = zeros(3,Li);
pos_mech = zeros(3,Li);
lla_mech = zeros(3,Li);

% initialization
pos(:,1) = lla2ecef(user_lla(1,:))';
vel(:,1) = gps.ecef_vel(:,1);
lla(:,1) = ecef2lla(pos(:,1)')';
pos_mech(:,1) = pos(:,1);
vel_mech(:,1) = vel(:,1);
lla_mech(:,1) = lla(:,1);

% NAV to ECEF transformation (Groves 2.150)
Rot_B_N = utils.euler2rot(psi_init);
[Rot_N_E, Rot_E_N] = utils.ned2ecef([deg2rad(lla(1:2,1)); lla(3,1)]);
C(:,:,1) = Rot_N_E * utils.euler2rot(psi_init);
att(:,1) = psi_init;
C_mech(:,:,1) = C(:,:,1);
att_mech(:,1) = psi_init;

% initial gps position
svInUse = gps.svInUse;
psr = gps.psr(svInUse,1);
dopp = gps.dopp(svInUse,1);

gps_pos = zeros(3,Lg);
gps_vel = zeros(3,Lg);
gps_lla = zeros(3,Lg);
[gps_pos(:,1),~,~,~,~] = utils.gnssPos(gps.x_sv(svInUse,:,1), psr, 1, [0;0;0;0]);
[gps_vel(:,1), ~] = utils.gnssVel(gps.x_sv(svInUse,:,1), gps.v_sv(svInUse,:,1), gps_pos(:,1), dopp);
gps_lla(:,1) = ecef2lla(gps_pos(:,1)')';

% initial ekf correction
if strcmp(kalmanType, 'tight')
    x = zeros(17,Lg);
    P = zeros(17,17,Lg);
    x(:,1) = [zeros(9,1); b_ad; b_gd; 0; 0];
    P(:,:,1) = eye(17);
    R = diag([1.0.*ones(length(psr),1); 0.1.*ones(length(dopp),1)]);
    Q = diag([3e0.*ones(3,1); 8e4.*ones(3,1); 0;0;0; 1e0.*ones(3,1); 6e-1.*ones(3,1); 4e-2; 1e-3]);
    [x(:,1),P(:,:,1)] = ekfTight(imu.w_ib_b(:,1)-b_gs, imu.f_ib_b(:,1)-b_as, C(:,:,1), pos(:,1), vel(:,1), ...
                psr, dopp, gps.x_sv(svInUse,:,1), gps.v_sv(svInUse,:,1), ...
                x(:,1), P(:,:,1), Q, R, 1, 'meas');
elseif strcmp(kalmanType, 'loose')
    x = zeros(15,Lg);
    P = zeros(15,15,Lg);
    x(:,1) = [zeros(9,1); b_ad; b_gd];
    P(:,:,1) = eye(15);
    R = diag([1.0.*ones(3,1); 0.1.*ones(3,1)]);
    Q = diag([3e0.*ones(3,1); 2e4.*ones(3,1); 0;0;0; 1e0.*ones(3,1); 5e-1.*ones(3,1)]);
    [x(:,1),P(:,:,1)] = ekfLoose(imu.w_ib_b(:,1)-b_gs, imu.f_ib_b(:,1)-b_as, C(:,:,1), pos(:,1), vel(:,1), ...
                gps_pos(:,1), gps_vel(:,1), x(:,1), P(:,:,1), Q, R, 1, 'meas');
end


%% SIMULATION

for i = 2:Li
    % imu time
    dti = imu.time(i) - imu.time(i-1);

    % imu measurement correction
    wb = imu.w_ib_b(:,i) - b_gd - b_gs;
    fb = imu.f_ib_b(:,i) - b_ad - b_as;

    % imu mechanization
    [pos(:,i), vel(:,i), C(:,:,i)] = ...
        imuMechECEF(wb, fb, pos(:,i-1), vel(:,i-1), C(:,:,i-1), dti);

    % save mechanization
    [pos_mech(:,i), vel_mech(:,i), C_mech(:,:,i)] = ...
        imuMechECEF(wb, fb, pos_mech(:,i-1), vel_mech(:,i-1), C_mech(:,:,i-1), dti);

    % check for new gps measurement
    g = find (gps.time > (imu.time(i-1)) & gps.time <= (imu.time(i)));

    % kalman update
    if (~isempty(g))
        
        % gps position update
        svInUse = gps.svInUse;
        psr = gps.psr(svInUse,g);
        dopp = gps.dopp(svInUse,g);
        svPos = gps.x_sv(svInUse,:,g);
        svVel = gps.v_sv(svInUse,:,g);
    
        % gps position and velocity
        gps_pos(:,g) = utils.gnssPos(svPos, psr, 1, [0;0;0;0]);
        gps_vel(:,g) = utils.gnssVel(svPos, svVel, gps_pos(:,g), dopp);
        gps_lla(:,g) = ecef2lla(gps_pos(:,g)')';

        % error state kalman filter
        if runKalman
            dtg = gps.gpsTime(g) - gps.gpsTime(g-1);
    
            % tightly coupled extended kalman filter
            if strcmp(kalmanType, 'tight')
                R = diag([1.0.*ones(length(psr),1); 0.01.*ones(length(dopp),1)]);
                [x(:,g),P(:,:,g)] = ekfTight(wb, fb, C(:,:,i), pos(:,i), vel(:,i), ...
                                             psr, dopp, svPos, svVel, ...
                                             x(:,g-1), P(:,:,g-1), Q, R, dtg, 'all');
            elseif strcmp(kalmanType, 'loose')
                R = diag([1.0.*ones(3,1); 0.1.*ones(3,1)]);
                [x(:,g),P(:,:,g)] = ekfLoose(wb, fb, C(:,:,i), pos(:,i), vel(:,i), ...
                            gps_pos(:,g), gps_vel(:,g), x(:,g-1), P(:,:,g-1), Q, R, dtg, 'all');
            end
    
            % Attitude correction (Groves 14.7)
            C(:,:,i) = (eye(3) - utils.vec2skew(x(1:3,g))) * C(:,:,i);
    
            % Velocity correction (Groves 14.8)
            vel(:,i)  = vel(:,i) - x(4:6,g);
    
            % Position correction (Groves 14.9)
            pos(:,i)  = pos(:,i) - x(7:9,g);
    
            % Biases estimation
            b_ad = b_ad + x(10:12,g);
            b_gd = b_ad + x(13:15,g);
    
            g = [];
        end

    end

    % attitude output in NED (euler angles)
    lla(:,i) = ecef2lla(pos(:,i)')';
    [~, Rot_E_N] = utils.ned2ecef([deg2rad(lla(1:2,i)); lla(3,i)]);
    att(:,i) = utils.rot2euler(Rot_E_N * C(:,:,i));
    lla_mech(:,i) = ecef2lla(pos_mech(:,i)')';
    [~, Rot_E_N] = utils.ned2ecef([deg2rad(lla_mech(1:2,i)); lla_mech(3,i)]);
    att_mech(:,i) = utils.rot2euler(Rot_E_N * C_mech(:,:,i));
end


%% PLOTTING

% f = figure(Units='normalized', Position=[1.1, 0.26, 1.5, 0.5]);
f = figure(Units='normalized', Position=[2.9, 0.375, 1.5, 0.49]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Parent=tbs, Title='Geoplot');
tab(2) = uitab(Parent=tbs, Title='3D Plot');
tab(3) = uitab(Parent=tbs, Title='Att. Error');
tab(4) = uitab(Parent=tbs, Title='Vel. Error');
tab(5) = uitab(Parent=tbs, Title='Vel. Err Zoom');
tab(6) = uitab(Parent=tbs, Title='Pos. Error');
tab(7) = uitab(Parent=tbs, Title='Pos. Err Zoom');
tab(8) = uitab(Parent=tbs, Title='Bias');

% GEOPLOT
geoaxes(Parent=tab(1));
geoplot(user_lla(:,1), user_lla(:,2), 'b', LineWidth=2.5);
hold on;
geoplot(gps_lla(1,:), gps_lla(2,:), 'g.', MarkerSize=10, LineWidth=5);
geoplot(lla(1,:), lla(2,:), 'r--', LineWidth=2);
geobasemap streets;
legend("REF", "GPS", "GNSS-INS", Location="northoutside", Orientation='horizontal');

% 3D PLOT
axes(Parent=tab(2));
hold on;
plot3(user_xyz(:,1), user_xyz(:,2), user_xyz(:,3), 'b', LineWidth=2.5);
plot3(gps_pos(:,1), gps_pos(:,2), gps_pos(:,3), 'g.', MarkerSize=10, LineWidth=5);
plot3(pos(:,1), pos(:,2), pos(:,3), 'r--', LineWidth=2);
grid on;
view(30, -30);
xlabel('x')
ylabel('y')
zlabel('z')

% ATTITUDE ERROR
axes(Parent=tab(3));
sgtitle("attitude");
subplot(3,1,1);
hold on;
plot(imu.time, rad2deg(ref.psi_n2b(3,:) - ref.psi_n2b(3,:)), 'b', LineWidth=2.5);
plot(imu.time, wrapTo180(rad2deg(att(3,:) - ref.psi_n2b(3,:))), 'r--', LineWidth=2);
plot(imu.time, wrapTo180(rad2deg(att_mech(3,:) - ref.psi_n2b(3,:))), 'k.', LineWidth=2);
% plot(t(1:100:10001), wrapTo180(rad2deg(x(3,:))), 'c--', LineWidth=2);
% legend("REF", "GNSS-INS", "IMU", "KF", Location="northoutside", Orientation='horizontal');
legend("REF", "GNSS-INS", "IMU", Location="northoutside", Orientation='horizontal');

subplot(3,1,2);
hold on;
plot(imu.time, rad2deg(ref.psi_n2b(2,:) - ref.psi_n2b(2,:)), 'b', LineWidth=2.5);
plot(imu.time, rad2deg(att(2,:) - ref.psi_n2b(2,:)), 'r--', LineWidth=2);
plot(imu.time, rad2deg(att_mech(2,:) - ref.psi_n2b(2,:)), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(2,:)), 'c--', LineWidth=2);

subplot(3,1,3);
hold on;
plot(imu.time, rad2deg(ref.psi_n2b(1,:) - ref.psi_n2b(1,:)), 'b', LineWidth=2.5);
plot(imu.time, rad2deg(att(1,:) - ref.psi_n2b(1,:)), 'r--', LineWidth=2);
plot(imu.time, rad2deg(att_mech(1,:) - ref.psi_n2b(1,:)), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(3,:)), 'c--', LineWidth=2);

% VELOCITY ERROR
axes(Parent=tab(4));
sgtitle("velocity");
subplot(3,1,1);
hold on;
plot(t, user_vel(1,:)-user_vel(1,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(1,:)-user_vel(1,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(1,:)-user_vel(1,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(1,:)-user_vel(1,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(4,:)), 'c--', LineWidth=2);
% legend("REF", "GPS", "GNSS-INS", "IMU", "KF", Location="northoutside", Orientation='horizontal');
legend("REF", "GPS", "GNSS-INS", "IMU", Location="northoutside", Orientation='horizontal');

subplot(3,1,2);
hold on;
plot(t, user_vel(2,:)-user_vel(2,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(2,:)-user_vel(2,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(2,:)-user_vel(2,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(2,:)-user_vel(2,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(5,:)), 'c--', LineWidth=2);

subplot(3,1,3);
hold on;
plot(t, user_vel(3,:)-user_vel(3,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(3,:)-user_vel(3,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(3,:)-user_vel(3,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(3,:)-user_vel(3,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(6,:)), 'c--', LineWidth=2);

% VELOCITY ERROR ZOOM
axes(Parent=tab(5));
sgtitle("velocity");
subplot(3,1,1);
hold on;
plot(t, user_vel(1,:)-user_vel(1,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(1,:)-user_vel(1,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(1,:)-user_vel(1,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(1,:)-user_vel(1,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(4,:)), 'c--', LineWidth=2);
% legend("REF", "GPS", "GNSS-INS", "IMU", "KF", Location="northoutside", Orientation='horizontal');
legend("REF", "GPS", "GNSS-INS", "IMU", Location="northoutside", Orientation='horizontal');
ylim([-0.6,0.6]);

subplot(3,1,2);
hold on;
plot(t, user_vel(2,:)-user_vel(2,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(2,:)-user_vel(2,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(2,:)-user_vel(2,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(2,:)-user_vel(2,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(5,:)), 'c--', LineWidth=2);
ylim([-0.6,0.6]);

subplot(3,1,3);
hold on;
plot(t, user_vel(3,:)-user_vel(3,:), 'b', LineWidth=2.5);
plot(gps.time, gps_vel(3,:)-user_vel(3,1:100:10001), 'g.', MarkerSize=10, LineWidth=5);
plot(imu.time, vel(3,:)-user_vel(3,:), 'r--', LineWidth=2);
plot(imu.time, vel_mech(3,:)-user_vel(3,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(6,:)), 'c--', LineWidth=2);
ylim([-0.6,0.6]);

% POSITION ERROR
axes(Parent=tab(6));
sgtitle("position error");
subplot(3,1,1);
hold on;
plot(t, user_xyz(1,:)-user_xyz(1,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(1,:)-user_xyz(1,1:100:10001), 'g*');
plot(imu.time, pos(1,:)-user_xyz(1,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(1,:)-user_xyz(1,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(7,:)), 'c--', LineWidth=2);
% legend("REF", "GPS", "GNSS-INS", "IMU", "KF", Location="northoutside", Orientation='horizontal');
legend("REF", "GPS", "GNSS-INS", "IMU", Location="northoutside", Orientation='horizontal');

subplot(3,1,2);
hold on;
plot(t, user_xyz(2,:)-user_xyz(2,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(2,:)-user_xyz(2,1:100:10001), 'g*');
plot(imu.time, pos(2,:)-user_xyz(2,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(2,:)-user_xyz(2,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(8,:)), 'c--', LineWidth=2);

subplot(3,1,3);
hold on;
plot(t, user_xyz(3,:)-user_xyz(3,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(3,:)-user_xyz(3,1:100:10001), 'g*');
plot(imu.time, pos(3,:)-user_xyz(3,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(3,:)-user_xyz(3,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(9,:)), 'c--', LineWidth=2);

% POSITION ERROR ZOOM
axes(Parent=tab(7));
sgtitle("position error");
subplot(3,1,1);
hold on;
plot(t, user_xyz(1,:)-user_xyz(1,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(1,:)-user_xyz(1,1:100:10001), 'g*');
plot(imu.time, pos(1,:)-user_xyz(1,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(1,:)-user_xyz(1,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(7,:)), 'c--', LineWidth=2);
% legend("REF", "GPS", "GNSS-INS", "IMU", "KF", Location="northoutside", Orientation='horizontal');
legend("REF", "GPS", "GNSS-INS", "IMU", Location="northoutside", Orientation='horizontal');
ylim([-3,3]);

subplot(3,1,2);
hold on;
plot(t, user_xyz(2,:)-user_xyz(2,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(2,:)-user_xyz(2,1:100:10001), 'g*');
plot(imu.time, pos(2,:)-user_xyz(2,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(2,:)-user_xyz(2,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(8,:)), 'c--', LineWidth=2);
ylim([-3,3]);

subplot(3,1,3);
hold on;
plot(t, user_xyz(3,:)-user_xyz(3,:), 'b', LineWidth=2.5);
plot(gps.time, gps_pos(3,:)-user_xyz(3,1:100:10001), 'g*');
plot(imu.time, pos(3,:)-user_xyz(3,:), 'r--', LineWidth=2);
plot(imu.time, pos_mech(3,:)-user_xyz(3,:), 'k.', LineWidth=2);
% plot(t(1:100:10001), rad2deg(x(9,:)), 'c--', LineWidth=2);
ylim([-3,3]);

% BIAS
axes(Parent=tab(8));
if strcmp(kalmanType, 'tight')
    subplot(2,2,1);
    plot(t(1:100:10001), x(10:12,:));
    title("Accelerometer Bias");
    
    subplot(2,2,2);
    plot(t(1:100:10001), x(13:15,:));
    title("Gyroscope Bias");
    
    subplot(2,2,3);
    plot(t(1:100:10001), x(16,:));
    title("Clock Bias");
    
    subplot(2,2,4);
    plot(t(1:100:10001), x(17,:));
    title("Clock Drift");
elseif strcmp(kalmanType, 'loose')
        subplot(2,1,1);
    plot(t(1:100:10001), x(10:12,:));
    title("Accelerometer Bias");
    
    subplot(2,1,2);
    plot(t(1:100:10001), x(13:15,:));
    title("Gyroscope Bias");
end


