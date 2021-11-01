function y = convolution(x, h)
    N = length(x);
    M = length(h);
    Ny = N + M -1;
    y = zeros(1,Ny);
    for i = 1:N
        for k = 1:M
            y(i+k-1) = y(i+k-1) + h(k)*x(i);
        end
%         y(i:i+M-1) = (y(i:i+M-1) - mean(y(i:i+M-1))) / std(y(i:i+M-1));
    end
end

