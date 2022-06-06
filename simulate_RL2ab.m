% Anne Collins 2022
% Code for Delevich et al 2022
% RL modeling
% this necessitates having Stan and MatlabStan installed



function [ttc]=simulate_RL2ab(RLData,Q0,par)
% returns simuluated data in the same format as input data RLData, given
% model parameters par, and initialization Q-values Q0
% simulated data is formatted as trials to criterion ttc.


epsilon=.00001;
for si=unique(RLData.Subject)
        alpha = par(si,1:2);
    beta=par(si,4:5);
    for iter=1:100
    Q=Q0;
    ph = 1;
    keepgoing=1;
    t=0;
    while keepgoing
        t=t+1;
        proba=epsilon/4+(1-epsilon)*exp(beta(ph)*Q)./(sum(exp(beta(ph)*Q)));
        cdf = [0 cumsum(proba)];
        a=find(cdf<rand);a=a(end);
        Reward(t)=(a==1);
        Phase(t)=ph;
        Q(a)=Q(a)+alpha(Phase(t))*(Reward(t)-Q(a));
        
        %[t proba Q Reward(t) ph]
        if Phase(t)==1 & t>=10
            criterion=sum(Reward(t-9:t));
            if criterion>=8
                ph =ph+1;
                ttc(si,1,iter)=t;
                criterion=0;
            end
        end
        if Phase(t)==2 & t>=ttc(si,1,iter)+10
            criterion=sum(Reward(t-9:t));
            if criterion>=8
                ph =ph+1;
                ttc(si,2,iter)=t-ttc(si,1,iter);
                keepgoing=0;
            end
        end
    end
    %boum
    end
end

ttc = mean(ttc,3);
end