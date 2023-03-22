function svData = svUnpack(struct, mode)

c = 299792458;
L = height(struct.measurements.L1.psr);
svData = cell(L,1);

for i = 1:L
    % determine usable satellites
    svInUse_L1 = ~isnan(struct.measurements.L1.psr(i,:));
    if strcmp(mode, 'sta')
        svInUse_L2 = ~isnan(struct.measurements.L2.psr(i,:));
    elseif strcmp(mode, 'dyn')
        svInUse_L2 = ~isnan(struct.measurements.L2.doppler(i,:));
    end

    % only parse if there are at least 4 SV
    if sum(svInUse_L1) > 3

        % pull in data
        svData{i}.svInUse = find(svInUse_L1 | svInUse_L2)';
        x = find(svInUse_L1 ~= svInUse_L2);   
        cond = svData{i}.svInUse==x;
        cond = sum(cond,2);
        svData{i}.svInUse = svData{i}.svInUse(~cond);

        svData{i}.gpsTime = struct.GPS_time.seconds(i);
        svData{i}.ephemeris = reshape(cell2mat(struct2cell(struct.ephem(svData{i}.svInUse))), [], length(svData{i}.svInUse))';
        svData{i}.L1_psr = struct.measurements.L1.psr(i,svData{i}.svInUse)';
        svData{i}.L1_dopp = struct.measurements.L1.doppler(i,svData{i}.svInUse)' .* -(c/1575.42e6); % Hz to m/s
        svData{i}.L1_car = struct.measurements.L1.carrier_phase(i,svData{i}.svInUse)';
        if strcmp(mode, 'sta')
            svData{i}.L2_psr = struct.measurements.L2.psr(i,svData{i}.svInUse)';
        end
        svData{i}.L2_dopp = struct.measurements.L2.doppler(i,svData{i}.svInUse)' .* -(c/1227.60e6); % Hz to m/s
        svData{i}.L2_car = struct.measurements.L2.carrier_phase(i,svData{i}.svInUse)';
  
        % create sv positions
        for j = 1:length(svData{i}.svInUse)
            transitTime = svData{i}.L1_psr(j) ./ c;
            transmitTime = svData{i}.gpsTime - transitTime;
            [svData{i}.svPos(j,:), svData{i}.svVel(j,:), svData{i}.clkCorr(j)] = ...
                utils.svPVTBevly(struct.ephem(svData{i}.svInUse(j)), transmitTime, transitTime);
        end
        svData{i}.clkCorr = svData{i}.clkCorr';

        % correct psuedorange
% %         if x(:,1) ~= x(:,2)
%         if sum(svInUse_L2) ~= sum(svInUse_L1)
%             cond = svData{i}.svInUse==x;
%         else
%             cond = zeros(size(svData{i}.clkCorr));
%         end
        svData{i}.L1_psr = svData{i}.L1_psr + svData{i}.clkCorr.*c;
        if strcmp(mode, 'sta')
            svData{i}.L2_psr = svData{i}.L2_psr + svData{i}.clkCorr.*c;
        end

%         svs = find(svInUse_L2);
%         for j = 1:length(svData{i}.svInUse)
%             transitTime = svData{i}.L2_psr(j) ./ c;
%             transmitTime = svData{i}.gpsTime - transitTime;
%             [svData{i}.L2_pos(j,:), svData{i}.L2_vel(j,:), svData{i}.L2_clkCorr(j)] = ...
%                 utils.svPVTBevly(struct.ephem(svs(j)), transmitTime, transitTime);
%         end
%         svData{i}.L2_clkCorr = svData{i}.L2_clkCorr';

    end
end

end

