function G2i = caGen(s1, s2)

% initialize to all ones
G1 = ones(1,10);
G2 = ones(1,10);
G2i = zeros(1,1023);

for i = 1:1023

    % create ca code
    G2i(i) = bitSum(G2(s1), G2(s2));
    G2i(i) = bitSum(G2i(i), G1(10));
    
    % update G1
    newBit = bitSum(G1(3), G1(10));
    G1 = [newBit, G1(1:9)];

    % update G2
    newBit = 0;
    for j = [2,3,6,8,9,10]
        newBit = bitSum(newBit, G2(j));
    end
    G2 = [newBit, G2(1:9)];

end
% fprintf("0x%c%c%c%c\n", dec2hex(bin2dec(num2str(G2i(1:16)))));
% G2i(G2i == 0) = -1;

% XOR function
function s = bitSum(a,b)
    s = a+b;
    if s > 1
        s = 0;
    end
end

end