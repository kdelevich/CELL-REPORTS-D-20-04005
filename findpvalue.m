% Anne Collins 2022
% Code for Delevich et al 2022
% RL modeling
% this necessitates having Stan and MatlabStan installed



function [p] = findpvalue(X)

% returns the [p, CI95] the p-value against 0 and the 95 credible interval
X = sort(X);
signchange = [1;sign(X(2:end).*X(1:end-1))];
p = find(signchange<=0)/length(X);
if isempty(p)
    p=0;
else
p=p(1);
end
if p>.5
    p=1-p;
end
p=2*p;


CI95(1) = X(round(0.025*length(X)));
CI95(2) = X(round(0.975*length(X)));

p = [p CI95];
end