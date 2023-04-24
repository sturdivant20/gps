function [R, shifts] = correlation(seq1, seq2, mode, shifts)
    if size(seq1) ~= size(seq2)
        fprintf("ERROR: 'correlation' -> sequneces not the same size!");
        return;
    end

    N = length(seq1);

    if nargin < 4
        shifts = -(N-1):(N-1);
    end
    if nargin < 3
        mode = 'standard';
    end

    R = zeros([length(shifts),1]);

    for k = 1:length(shifts)
        seq2_ = circshift(seq2, shifts(k));
        % seq2_ = [seq2(N-shifts(k)+1:N), seq2(1:N-shifts(k))];
        % R(k) = 1/N * sum(seq1.*seq2_);
        R(k) = seq1*seq2_';
    end

    if strcmp(mode,'normalized')
        R = R./N;
    end
end

