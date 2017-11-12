function [ w1, w2 ] = compute_weight(P,K_min,K_max)

if P(1) < K_min(1) 
    w1 = 0;
    w2 = 1;
    return;
end

if P(1) > K_max(1)
    w1 = 1;
    w2 = 0;
    return;
end

a = (P(1) - K_min(1)) * (K_max(1) - K_min(1));
b = (P(2) - K_min(2)) * (K_max(2) - K_min(2));
c = (K_max(1) - K_min(1))^2 + (K_max(2) - K_min(2))^2;

w1 = abs(a + b)/c;
w2 = 1 - w1;
end
