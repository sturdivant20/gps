%% Unpack L1 GPS Data | Daniel Sturdivant
clc; clear; close all;

% choose wheter to output new data file
saveFile = true;

% constants
load("+data/RCVR_D1_data_new2.mat");
load("+data/IMU_new2.mat");
c = 299792458;
lenG = length(RCVR_D1.GPS_time.seconds);
lenI = length(IMU.time.full_time);
L1 = 1575.42e6;

% gps data to retreive
week = RCVR_D1.GPS_time.week;
svInUse = boolean(zeros(30,lenG));
timeG = zeros(1,lenG);
psr = zeros(30,lenG);
dopp = zeros(30,lenG);
car = zeros(30,lenG);
cn0 = zeros(30,lenG);
x_sv = zeros(30,3,lenG);
v_sv = zeros(30,3,lenG);
lla = zeros(3,lenG);

% unpack gps data
for i = 1:lenG

    % determine usable satellites
    svInUse(:,i) = (~isnan(RCVR_D1.measurements.L1.psr(i,:)) ...
                 & ~isnan(RCVR_D1.measurements.L1.doppler(i,:)) ...
                 & ~isnan(RCVR_D1.measurements.L1.carrier_phase(i,:)))';
    k = find(svInUse(:,i));

    % pull in data
    timeG(i) = RCVR_D1.GPS_time.seconds(i) - RCVR_D1.GPS_time.seconds(1);
    psr(k,i) = RCVR_D1.measurements.L1.psr(i,svInUse(:,i))';
    dopp(k,i) = RCVR_D1.measurements.L1.doppler(i,svInUse(:,i))' .* -(c/L1);
    car(k,i) = RCVR_D1.measurements.L1.carrier_phase(i,svInUse(:,i))' .* (c/L1);
    cn0(k,i) = RCVR_D1.measurements.L1.carrier_to_noise(i,svInUse(:,i))';

    % create satellite positions and velocities
    clkCorr = zeros(sum(svInUse(:,i)),1);
    for j = 1:sum(svInUse(:,i))
        transitTime = psr(k(j),i) ./ c;
        transmitTime = RCVR_D1.GPS_time.seconds(i) - transitTime;
        [x_sv(k(j),:,i), v_sv(k(j),:,i), clkCorr(j)] = ...
            data.calc_gps_sv_pos(RCVR_D1.ephem(k(j)), transmitTime, transitTime);
    end

    % correct psuedorange
    psr(k,i) = psr(k,i) + clkCorr .* c;

    % true posistion
    lla(:,i) = [RCVR_D1.true_pos.lat(i); RCVR_D1.true_pos.lon(i); RCVR_D1.true_pos.alt(i)];

end


% imu data to retreive
timeI = zeros(1,lenI);
w_ib_b = zeros(3,lenI);
f_ib_b = zeros(3,lenI);
quat = zeros(4,lenI);
euler = zeros(3,lenI);

% unpack imu data
for i = 1:lenI

    % pull in data
    timeI(i) = IMU.time.full_time(i) - IMU.time.full_time(1);
    w_ib_b(:,i) = [IMU.angular_velocity.X(i); IMU.angular_velocity.Y(i); IMU.angular_velocity.Z(i)];
    f_ib_b(:,i) = [IMU.linear_acceleration.X(i); IMU.linear_acceleration.Y(i); IMU.linear_acceleration.Z(i)];
    quat(:,i) = [IMU.orientation.X(i); IMU.orientation.Y(i); IMU.orientation.Z(i); IMU.orientation.W(i)];
    euler(:,i) = utils.quat2euler(quat(:,i));

end


% create and save struct
if saveFile == true

    % gps data
    gps.svInUse = svInUse;
    gps.gpsTime = timeG;
    gps.psr = psr;
    gps.dopp = dopp;
    gps.car = car;
    gps.cn0 = cn0;
    gps.x_sv = x_sv;
    gps.v_sv = v_sv;
    gps.lla = lla;

    % imu data
    imu.time = timeI;
    imu.w_ib_b = w_ib_b;
    imu.f_ib_b = f_ib_b;
    imu.quat = quat;
    imu.euler = euler;

    save("+data/data.mat", "gps", "imu");
end