function svData = svUnpack(struct, time)

c = 299792458;
L = height(struct.measurements.L1.psr);
svData = cell(L,1);

for i = 1:L
    % determine usable satellites
    svInUse = ~isnan(struct.measurements.L1.psr(i,:));
    if(i == 121)
        i
    end

    % only parse if there are at least 4 SV
    if sum(svInUse) > 3

        % pull in data
        svData{i}.svInUse = find(svInUse)';
        svData{i}.gpsTime = struct.GPS_time.seconds(i);
        svData{i}.ephemeris = cell2mat(struct2cell(struct.ephem(svInUse)))';
        svData{i}.psr = struct.measurements.L1.psr(i,svInUse)';
        svData{i}.dopp = struct.measurements.L1.doppler(i,svInUse)' .* -(c/1575.42e6); % Hz to m/s
    
        if (nargin > 1) && (i >= 1068)
            svData{i}.gpsTime = svData{i}.gpsTime + 1;
        end
        % check satellite data
        k = 0;
        for j = 1:length(svData{i}.svInUse)
            if isnan(svData{i}.ephemeris(j-k,2))
                svData{i}.svInUse(j-k) = [];
                svData{i}.psr(j-k) = [];
                svData{i}.dopp(j-k) = [];
                svData{i}.ephemeris(j-k,:) = [];
                k = k + 1;
            end
        end
    
        % create sv positions
        for j = 1:length(svData{i}.svInUse)
            transitTime = svData{i}.psr(j) ./ c;
            transmitTime = svData{i}.gpsTime - transitTime;
            [svData{i}.pos(j,:), svData{i}.vel(j,:), svData{i}.clkCorr(j)] = ...
                svPVTBevly(struct.ephem(svData{i}.svInUse(j)), transmitTime, transitTime);
        end
        svData{i}.clkCorr = svData{i}.clkCorr';
    %     transitTime = svData{i}.psr ./ c;
    %     transmitTime = svData{i}.gpsTime - transitTime;
    %     [svData{i}.clkCorr, svData{i}.pos, svData{i}.vel] = ...
    %         svPVT(svData{i}.ephemeris, transmitTime, transitTime);
    
        % correct psuedorange
        svData{i}.psr_corr = svData{i}.psr + svData{i}.clkCorr.*c;
    
        % elevation mask
        lla = [32.58632572006568, -85.49423209499069, 250];
        [svData{i}.az, svData{i}.el, ~] = ecef2aer(svData{i}.pos(:,1), svData{i}.pos(:,2), svData{i}.pos(:,3), ...
                                lla(1), lla(2), lla(3), wgs84Ellipsoid("meter"));
    %     elMask = abs(svData{i}.el) >= 10;
    %     svData{i}.svInUse = svData{i}.svInUse(elMask);
    %     svData{i}.pos = svData{i}.pos(elMask, :);
    %     svData{i}.vel = svData{i}.vel(elMask, :);
    %     svData{i}.psr_corr = svData{i}.psr_corr(elMask);
    %     svData{i}.dopp = svData{i}.dopp(elMask);
    %     svData{i}.ephemeris = svData{i}.ephemeris(elMask, :);
    %     svData{i}.clkCorr = svData{i}.clkCorr(elMask);
    %     svData{i}.az = svData{i}.az(elMask);
    %     svData{i}.el = svData{i}.el(elMask);
    

        % create user positions and velocities
        [svData{i}.x,svData{i}.b,svData{i}.P,svData{i}.DOP,svData{i}.i] = ...
            gnssPos(svData{i}.pos, svData{i}.psr_corr, 0.5);
        [svData{i}.xDot, svData{i}.bDot] = ...
            gnssVel(svData{i}.pos, svData{i}.vel, svData{i}.x, svData{i}.dopp);
    
        % lla postition
        svData{i}.lla = ecef2lla(svData{i}.x');

    end
end

end

