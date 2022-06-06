% Anne Collins 2022
% Code for Delevich et al 2022
% RL modeling
% this necessitates having Stan and MatlabStan installed


function [ps,entropy,finalQs]=computeLLH_RL2ab(RLData,Q0,par)

epsilon=.00001;
k = 0;
for si=unique(RLData.Subject)
    Choice=RLData.Choice(RLData.Subject==si)';
    Reward=RLData.Reward(RLData.Subject==si)';
    Phase=RLData.Phase(RLData.Subject==si)';
    beta = par(si,1:2);
    alpha=par(si,3:4);
    Q=Q0;
    for t=1:length(Choice)
        k=k+1;
        if t>1 & Phase(t)>1 & Phase(t-1)==1
            proba=exp(beta(1)*Q)./(sum(exp(beta(1)*Q)));
            entropy(si,1) = sum(proba.*log(proba));
            proba=exp(beta(2)*Q)./(sum(exp(beta(2)*Q)));
            entropy(si,2)=sum(proba.*log(proba));
            finalQs(si,:)=Q;
        end
        ps(k) = epsilon/4+(1-epsilon)/(sum(exp(beta(Phase(t))*(Q-Q(Choice(t))))));
        Q(Choice(t))=Q(Choice(t))+alpha(Phase(t))*(Reward(t)-Q(Choice(t)));
    end
end
end