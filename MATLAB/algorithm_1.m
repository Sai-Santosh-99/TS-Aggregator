close all;
clear all;
clc;

prob = rand(1, 8);

horizon = 10000;
index = 1:1:horizon;

n_tests = 50;

regret_ts = zeros(n_tests, horizon+1);
reward_ts = zeros(n_tests, horizon+1);

for x=1:n_tests
    
    disp(x);
    
    X = ones(1,length(prob));
    T = ones(1,length(prob));
    Q = zeros(1, length(prob));
    N = length(prob);
    
    for i=1:horizon
        
        for j=1:length(prob)
            Q(j) = machine(X(j), T(j));
        end

        m = max(Q);
        I = find(Q == m);
        r = rand(1, length(prob));

        if r(I(1)) < prob(I(1))
            X(I(1)) = X(I(1)) + 1;
            reward_ts(x, i+1) = reward_ts(x, i) + 1;
        else
            reward_ts(x, i+1) = reward_ts(x, i);
        end

        regret_ts(x, i+1) = regret_ts(x, i) + max(prob) - prob(I(1));

        T(I(1)) = T(I(1)) + 1;
        N = N + 1;
    end
end

function Q = machine(X, T) 
    val = rand(1,T);
    val = sort(val);
    Q = val(X);
end
