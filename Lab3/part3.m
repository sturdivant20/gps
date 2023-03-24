%% GPS LAB 3 Part 3 | Daniel Sturdivant & Andrew Weir
clc; clear; close all;
fprintf("<strong>PART 3\n</strong>");

c = 299792458;  % speed of light
L1 = 1575.42e6;
L2 = 1227.60e6;
lambda1 = c / L1;
lambda2 = c / L2;

load("RCVR_D1.mat");
r1.L = length(D1);

r1.gpsTime = zeros(r1.L,1);
r1.L1psr = nan(32,r1.L);
r1.L1psr_var = nan(32,r1.L);
r1.L1dopp = nan(32,r1.L);
r1.L1car = nan(32,r1.L);
r1.L1car_var = nan(32,r1.L);
r1.L2psr = nan(32,r1.L);
r1.L2psr_var = nan(32,r1.L);
r1.L2dopp = nan(32,r1.L);
r1.L2car = nan(32,r1.L);
r1.L2car_var = nan(32,r1.L);

for i = 1:r1.L
    r1.gpsTime(i) = D1{i}.gpsTime;
    svTest = D1{i}.svInUse;

    r1.L1psr(svTest,i) = D1{i}.L1_psr;
    r1.L1psr_var(svTest,i) = D1{i}.L1_psr_var;
    r1.L1dopp(svTest,i) = D1{i}.L1_dopp;
    r1.L1car(svTest,i) = D1{i}.L1_car;
    r1.L1car_var(svTest,i) = D1{i}.L1_car_var;

    r1.L2psr(svTest,i) = D1{i}.L2_psr;
    r1.L2psr_var(svTest,i) = D1{i}.L2_psr_var;
    r1.L2dopp(svTest,i) = D1{i}.L2_dopp;
    r1.L2car(svTest,i) = D1{i}.L2_car;
    r1.L2car_var(svTest,i) = D1{i}.L2_car_var;
end

load("RCVR_D2.mat");
r2.L = length(D2);

r2.gpsTime = zeros(r2.L,1);
r2.L1psr = nan(32,r2.L);
r2.L1psr_var = nan(32,r2.L);
r2.L1dopp = nan(32,r2.L);
r2.L1car = nan(32,r2.L);
r2.L1car_var = nan(32,r1.L);
r2.L2psr = nan(32,r2.L);
r2.L2psr_var = nan(32,r2.L);
r2.L2dopp = nan(32,r2.L);
r2.L2car = nan(32,r1.L);
r2.L2car_var = nan(32,r1.L);

for i = 1:r2.L
    r2.gpsTime(i) = D2{i}.gpsTime;
    svTest = D2{i}.svInUse;

    r2.L1psr(svTest,i) = D2{i}.L1_psr;
    r2.L1psr_var(svTest,i) = D2{i}.L1_psr_var;
    r2.L1dopp(svTest,i) = D2{i}.L1_dopp;
    r2.L1car(svTest,i) = D2{i}.L1_car;
    r2.L1car_var(svTest,i) = D2{i}.L1_car_var;

    r2.L2psr(svTest,i) = D2{i}.L2_psr;
    r2.L2psr_var(svTest,i) = D2{i}.L2_psr_var;
    r2.L2dopp(svTest,i) = D2{i}.L2_dopp;
    r2.L2car(svTest,i) = D2{i}.L2_car;
    r2.L2car_var(svTest,i) = D2{i}.L2_car_var;
end


% only use portions where times are the same
r1.Start = find(r1.gpsTime == max([r1.gpsTime(1), r2.gpsTime(1)]));
r2.Start = find(r2.gpsTime == max([r1.gpsTime(1), r2.gpsTime(1)]));
r1.Finish = find(r1.gpsTime == min([r1.gpsTime(end), r2.gpsTime(end)]));
r2.Finish = find(r2.gpsTime == min([r1.gpsTime(end), r2.gpsTime(end)]));

r1.L = r1.Finish - r1.Start + 1;
r1.gpsTime = r1.gpsTime(r1.Start:r1.Finish);
r1.L1psr = r1.L1psr(:, r1.Start:r1.Finish);
r1.L1psr_var = r1.L1psr_var(:, r1.Start:r1.Finish);
r1.L1dopp = r1.L1dopp(:, r1.Start:r1.Finish);
r1.L1car = r1.L1car(:, r1.Start:r1.Finish);
r1.L2psr = r1.L2psr(:, r1.Start:r1.Finish);
r1.L2psr_var = r1.L2psr_var(:, r1.Start:r1.Finish);
r1.L2dopp = r1.L2dopp(:, r1.Start:r1.Finish);
r1.L2car = r1.L2car(:, r1.Start:r1.Finish);
r1.L2car_var = r1.L2car_var(:, r1.Start:r1.Finish);

r2.L = r2.Finish - r2.Start + 1;
r2.gpsTime = r2.gpsTime(r2.Start:r2.Finish);
r2.L1psr = r2.L1psr(:, r2.Start:r2.Finish);
r2.L1psr_var = r2.L1psr_var(:, r2.Start:r2.Finish);
r2.L1dopp = r2.L1dopp(:, r2.Start:r2.Finish);
r2.L1car = r2.L1car(:, r2.Start:r2.Finish);
r2.L2psr = r2.L2psr(:, r2.Start:r2.Finish);
r2.L2psr_var = r2.L2psr_var(:, r2.Start:r2.Finish);
r2.L2dopp = r2.L2dopp(:, r2.Start:r2.Finish);
r2.L2car = r2.L2car(:, r2.Start:r2.Finish);
r2.L2car_var = r2.L2car_var(:, r2.Start:r2.Finish);


%% Part A
fprintf("<strong>\n(a)\n</strong>");

% position comparison
a_x1 = zeros(3,r1.L);
a_lla1 = zeros(3,r1.L);
a_x2 = zeros(3,r2.L);
a_lla2 = zeros(3,r2.L);
a_dist = zeros(1,r1.L);

k1 = false(1,r1.L);
k2 = false(1,r2.L);
for i = 1:r1.L
    svTest1 = D1{i+r1.Start-1}.svInUse;
    if length(svTest1) > 3
        svPos1 = D1{i+r1.Start-1}.svPos;
        psr1 = r1.L1psr(svTest1, i);
        [a_x1(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos1, psr1, 0.5);
        a_lla1(:,i) = ecef2lla(a_x1(:,i)')';
    else
        k1(i) = 1;
    end

    svTest2 = D2{i+r2.Start-1}.svInUse;
    if length(svTest2) > 3
        svPos2 = D2{i+r2.Start-1}.svPos;
        psr2 = r2.L1psr(svTest2, i);
        [a_x2(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos2, psr2, 0.5);
        a_lla2(:,i) = ecef2lla(a_x2(:,i)')';
    else
        k2(i) = 1;
    end

    a_dist(i) = norm(a_x1(:,i) - a_x2(:,i));
end

a_x1(:,k1) = [];
a_x2(:,k2) = [];
a_lla1(:,k1) = [];
a_lla2(:,k2) = [];
a_dist(:,k1|k2) = [];

a_mu = mean(a_dist);
a_std = std(a_dist);

figure();
geoplot(a_lla1(1,:), a_lla1(2,:), 'o', LineWidth=2.5);
hold on;
geoplot(a_lla2(1,:), a_lla2(2,:), 'x', LineWidth=2.5);
legend("RCVR1", "RCVR2");
geobasemap satellite;
title("Regular Positioning")

fprintf("std = %f\n", a_std);
fprintf("mean = %f\n", a_mu);


%% Part B
fprintf("<strong>\n(b)\n</strong>");

% dgps solution
b_x = zeros(3,r2.L);
b_lla = zeros(3,r2.L);
b_dist = zeros(1,r1.L);

XX = zeros(4,r1.L);

k = false(1,r1.L);
sk = 0;
for i = 1:r1.L
    svTest1 = D1{i+r1.Start-1}.svInUse;
    svTest2 = D2{i+r2.Start-1}.svInUse;
    if sum(ismember(svTest1, svTest2)) > 3
        psr1 = r1.L1psr(svTest1(ismember(svTest1, svTest2)), i);
        psr2 = r2.L1psr(svTest2(ismember(svTest2, svTest1)), i);
        svPos = D2{i+r2.Start-1}.svPos(ismember(svTest2, svTest1),:); % 2 is the magic number

        svVel = D2{i+r2.Start-1}.svVel(ismember(svTest2, svTest1),:);
        dopp2 = r2.L1dopp(svTest2(ismember(svTest2, svTest1)), i);

        [b_0(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos, psr2, 0.5);
        [v_0(:,i), ~] = utils.gnssVel(svPos, svVel, b_0(:,i), dopp2);

        u = svPos - b_0(:,i)';
        r = sqrt(sum(u.^2,2));

        H = [u./r, ones(size(psr1))];
        y = psr1 - psr2;
        X = pinv(H)*y;

        if i == 1
            Q = H' * H;
            XX(:,1) = (Q^-1 * H'  * y);
        else
            Q = Q + (H'*H);
            dx = Q^-1 * H' * (y - H*XX(:,i-1-sk));
            XX(:,i) = XX(:,i-1-sk) + dx;
            sk = 0;
        end

        b_x(:,i) = b_0(:,i) + X(1:3);
        b_lla(:,i) = ecef2lla(b_x(:,i)')';
        b_dist(i) = norm(b_x(:,i) - b_0(:,i));
        b_x2(:,i) = b_0(:,i) + XX(1:3,i);
        b_lla2(:,i) = ecef2lla(b_x2(:,i)')';
        b_dist2(i) = norm(b_x2(:,i) - b_0(:,i));
    else
        k(i) = 1;
        sk = sum(k);
    end
end

b_0(:,k) = [];
v_0(:,k) = [];

b_x(:,k) = [];
b_lla(:,k) = [];
b_dist(:,k) = [];
b_x2(:,k) = [];
b_lla2(:,k) = [];
b_dist2(:,k) = [];

b_mu = mean(b_dist);
b_std = std(b_dist);
b_mu2 = mean(b_dist2);
b_std2 = std(b_dist2);

figure();
geoplot(a_lla2(1,:), a_lla2(2,:), 'o', LineWidth=2.5);
hold on;
geoplot(b_lla(1,:), b_lla(2,:), 'x', LineWidth=2.5);
geoplot(b_lla2(1,:), b_lla2(2,:), '.', LineWidth=2.5);
legend("RCVR2", "RCVR1 DGPS", "RCVR1 Recursive DGPS");
geobasemap satellite;
title("DGPS Positioning")

fprintf("std = %f\n", b_std);
fprintf("mean = %f\n", b_mu);
fprintf("std2 = %f\n", b_std2);
fprintf("mean2 = %f\n", b_mu2);


%% Part C
fprintf("<strong>\n(c)\n</strong>");

% carrier smoothed dgps
k = ones(32,1);
kk = false(1,r1.L);
M = [1/60,8,15].*60;

c_x = zeros(3,r2.L, 3);
c_lla = zeros(3,r2.L, 3);
c_dist = zeros(3,r1.L);

dP = zeros(32,r1.L);
dC = zeros(32,r1.L);

for m = 1:3
    for i = 1:r1.L
        svTest1 = D1{i+r1.Start-1}.svInUse;
        svTest2 = D2{i+r2.Start-1}.svInUse;
        svTest = svTest2(ismember(svTest2, svTest1));

        if sum(ismember(svTest1, svTest2)) > 3
            for j = 1:32
    
                if any(ismember(svTest,j))
        
                    dC(j,i) = lambda1*abs(r2.L1car(j,i)) - lambda1*abs(r1.L1car(j,i));
                    psr1 = r1.L1psr(j, i);
                    psr2 = r2.L1psr(j, i);
                    if (i == 1) || (k(j) == 1)
                       dP(j,i) = psr2 - psr1;
                    elseif k(j) < M(m)
                        ddC = dC(j,i) - dC(j,i-1);
                        if abs(ddC) > 10
                            ddC = 0;
                        end
                        dP(j,i) = (1/i)*(psr2 - psr1) + ((i-1)/i)*(dP(j,i-1) + ddC);
                    else 
                        ddC = dC(j,i) - dC(j,i-1);
                        if abs(ddC) > 10
                            ddC = 0;
                        end
                        dP(j,i) = (1/M(m))*(psr2 - psr1) + ((M(m)-1)/M(m))*(dP(j,i-1) + ddC);
                    end
                else
                    k(j) = 1;
                end

                k(j) = k(j) + 1;
        
            end

%             if i == 193
%                 pause
%             end

            % position
            svPos = D2{i+r2.Start-1}.svPos(ismember(svTest2, svTest1),:);
            psr = r2.L1psr(svTest, i);
            [tmp, ~, ~, ~, ~] = utils.gnssPos(svPos, psr, 0.5);
            u = svPos - tmp';
            r = sqrt(sum(u.^2,2));
            H = [u./r, ones(size(dP(svTest,i)))];
            y = dP(svTest,i);
            X = pinv(H)*y;
    
            c_x(:,i,m) = tmp + X(1:3);
            c_lla(:,i,m) = ecef2lla(c_x(:,i,m)')';
            c_dist(m,i) = norm(c_x(:,i,m) - tmp);

        else
            kk(i) = 1;
        end

    end
end

c_x(:,kk,:) = [];
c_lla(:,kk,:) = [];
c_dist(:,kk,:) = [];

c_mu1 = mean(c_dist(1,:));
c_std1 = std(c_dist(1,:));
c_mu2 = mean(c_dist(2,:));
c_std2 = std(c_dist(2,:));
c_mu3 = mean(c_dist(3,:));
c_std3 = std(c_dist(3,:));

figure();
geoplot(a_lla2(1,:), a_lla2(2,:), 'o', LineWidth=2.5);
hold on;
geoplot(c_lla(1,:,1), c_lla(2,:,1), 'x', LineWidth=3.5);
geoplot(c_lla(1,:,2), c_lla(2,:,2), 'x', LineWidth=1.5);
geoplot(c_lla(1,:,3), c_lla(2,:,3), 'x', LineWidth=0.5);
legend("RCVR2", "2 MIN", "8 MIN", "15 MIN");
geobasemap satellite;
title("Carrier Smoothed DGPS Positioning")

fprintf("std 2 min = %f\n", c_std1);
fprintf("mean 2 min = %f\n", c_mu1);
fprintf("std 8 min = %f\n", c_std2);
fprintf("mean 8 min = %f\n", c_mu2);
fprintf("std 15 min = %f\n", c_std3);
fprintf("mean 15 min = %f\n", c_mu3);


%% Part D
fprintf("<strong>\n(d)\n</strong>");

% rtk dgps (carrier)
d_x = zeros(3,r1.L);
d_lla = zeros(3,r1.L);
d_dist = zeros(1,r1.L);
% d_N = zeros(2*T,r1.L);
% d_N2 = zeros(2*T,r1.L);
% d_N_cov = zeros(2*T,2*T,r1.L);
kk = false(1,r1.L);
k = 1;

for i = 1:r1.L
    svTest1 = D1{i+r1.Start-1}.svInUse;
    svTest2 = D2{i+r2.Start-1}.svInUse;
    svTest = svTest2(ismember(svTest2, svTest1));

    if sum(ismember(svTest1, svTest2)) > 3
        T = length(svTest);
    
        svPos = D2{i+r2.Start-1}.svPos(ismember(svTest2, svTest1),:);
        psr = r2.L1psr(svTest, i);
        [tmp, ~, ~, ~, ~] = utils.gnssPos(svPos, psr, 0.5);
    
        u = svPos - tmp';
        r = sqrt(sum(u.^2,2));
        
        % recursive least squares for N
%         R = diag([r1.L1psr_var(svTest,i) + r2.L1psr_var(svTest,i); ...
%                   r1.L1car_var(svTest,i) + r2.L1car_var(svTest,i); ...
%                   r1.L2car_var(svTest,i) + r2.L2car_var(svTest,i)]);
    
        m = [r2.L1psr(svTest,i) - r1.L1psr(svTest,i); ...
             abs(r2.L1car(svTest,i).*lambda1) - abs(r1.L1car(svTest,i).*lambda1); ...
             abs(r2.L2car(svTest,i).*lambda2) - abs(r1.L2car(svTest,i).*lambda2)];
        
        G = [u./r, ones(T,1); ...
             u./r, ones(T,1); ...
             u./r, ones(T,1)];
    
        LAM = [zeros(T,T)              , zeros(T,T); ...
               lambda1.*diag(ones(T,1)), zeros(T,T); ...
               zeros(T,T)              , lambda2.*diag(ones(T,1))];
    
        L = null(G')';
    
%         RR = L*R^-1*L';
        HH = L*LAM;
        yy = L*m;
        d_N = (HH'*HH)^-1*HH'*yy;
        d_N_cov = (HH'*HH)^-1;
    
%         if (i == 1) || (k == 1)
%             Q = HH' * HH;
%             d_N = (Q^-1 * HH' * yy);
%             d_N_cov = eye(length(d_N)).*10;
%             k = k + 1;
% 
%         elseif isequal(svTest, svTestPrev)
%             Q = Q + (HH' * HH);
%             dx = Q^-1 * HH' * (yy - HH*d_N);
%             d_N = d_N2(:,1) + dx;
%             d_N_cov = Q^-1;
% 
% %         elseif any(~ismember(svTest, svTestPrev))
%         else
%             Q = HH' * HH;
%             d_N = (Q^-1 * HH' * yy);
%             d_N_cov = Q^-1;
% 
% %         elseif length(svTest) < length(svTestPrev)
% %         else
% %             cond = find(repmat(~ismember(svTestPrev, svTest),2,1));
% %             Q(cond,:) = [];
% %             Q(:,cond) = [];
% %             d_N2(cond,:) = [];
% % 
% %             Q = Q + (HH' * HH);
% %             dx = Q^-1 * HH' * (yy - HH*d_N2(:,1));
% %             d_N = d_N2(:,1) + dx;
% %             d_N_cov = Q^-1;
% 
%         end
    
        % lambda method
        [d_N2,~,~,~,~,~,~] = LAMBDA(d_N, d_N_cov);
    
        % find position
        y = abs(r2.L1car(svTest,i)*lambda1) - abs(r1.L1car(svTest,i)*lambda1) - (lambda1.*round(d_N2(1:T,1)));
        H = [u./r, ones(T,1)];
        dx = pinv(H)*y;
        d_x(:,i) = tmp + dx(1:3);
        d_lla(:,i) = ecef2lla(d_x(:,i)')';
    
        d_dist(i) = norm(tmp - d_x(:,i));

        svTestPrev = svTest;
    else
        kk(i) = 1;
    end

end

d_x(:,kk) = [];
d_lla(:,kk) = [];
d_dist(:,kk) = [];
d_mu = mean(d_dist);
d_std = std(d_dist);

figure();
geoplot(a_lla2(1,:), a_lla2(2,:), 'o');
hold on;
geoplot(d_lla(1,:), d_lla(2,:), 'x');
legend("RCVR2", "RCVR1 RTK DGPS");
geobasemap satellite;
title("RTK DGPS Positioning")

fprintf("std rtk = %f\n", d_std);
fprintf("mean rtk = %f\n", d_mu);


%% Part E
fprintf("<strong>\n(e)\n</strong>");

% attitude check
lla0 = ecef2lla(b_0')';
[E1, N1, U1] = ecef2enu(b_0(1,:), b_0(2,:), b_0(3,:), lla0(1,:), lla0(2,:), lla0(3,:), wgs84Ellipsoid("meter"));
[E2, N2, U2] = ecef2enu(b_x2(1,:), b_x2(2,:), b_x2(3,:), lla0(1,:), lla0(2,:), lla0(3,:), wgs84Ellipsoid("meter"));
[VE1, VN1, VU1] = ecef2enuv(v_0(1,:), v_0(2,:), v_0(3,:), lla0(1,:), lla0(2,:));
L = norm(b_0 - b_x2);

speed = vecnorm([VE1;VN1;VU1],2,1);
VE1(speed > 30) = NaN;
VN1(speed > 30) = NaN;
VU1(speed > 30) = NaN;
VE1 = fillgaps(VE1, 10, 1);
VN1 = fillgaps(VN1, 10, 1);
VU1 = fillgaps(VU1, 10, 1);

course = wrapTo360(atan2d(VE1 , VN1));
% course(speed < 2) = 0;
psi = wrapTo360(atand((E2-E1) ./ (N2-N1))-15);

figure();
hold on;
plot(course,'.b', LineWidth=2);
plot(psi,'r', LineWidth=2);
title("Heading Comparison");
xlabel("Time [s]");
ylabel("Angle [deg]");
legend("Course", "Azimuth");


%% dist plot
figure();
hold on;
plot(a_dist, 'o', LineWidth=2.5);
plot(b_dist, 'o', LineWidth=1.5);
plot(c_dist(1,:), 'x', LineWidth=2);
plot(c_dist(2,:), 'x', LineWidth=2);
plot(c_dist(3,:), 'x', LineWidth=2);
plot(d_dist(1,:), 'square', LineWidth=2);
plot(b_dist2, '^', LineWidth=1);
legend("A", "B", "C 2 MIN", "C 8 MIN", "C 15 MIN", "D", "B Recursive")
title("Dynamic Receiver Positioning Difference")
xlabel("Time [s]");
ylabel("Position Difference [m]");
