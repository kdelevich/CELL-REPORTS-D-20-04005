% Anne Collins 2022
% Code for Delevich et al 2022
% RL modeling
% this necessitates having Stan and MatlabStan installed


function RL_data=FormatData
% returns the data in a structure appropriate for matlabstan modeling
%clear all

folder = '4C_10.29.21/';
names = {'RLSplitAlphaBeta_parameters_10.29.21_D1hm4di',
     'RLSplitAlphaBeta_parameters_10.29.21_D1mCherry',
     'RLSplitAlphaBeta_parameters_10.29.21_D2hm3dq',
     'RLSplitAlphaBeta_parameters_10.29.21_D2hm4di',
     'RLSplitAlphaBeta_parameters_10.29.21_D2mCherry',
     'RLSplitAlphaBeta_parameters_10.29.21_Panhm4di'};
% folder='Data/';
% names = {'d1hm4di_trialhistory','d2hm4di_trialhistory',...
%     'd2hm3dq_trialhistory','mcherry_trialhistory'};
%names = {'d1hm4di','mcherry'};
% 
% folder='Data/final4choicecohortdatawithmcherrysplitbyd1crean/';
% names = {'RLSplitAlphaForget_parameters_D1hm4di',...
%     'RLSplitAlphaForget_parameters_D2hm4di',...
%     'RLSplitAlphaForget_parameters_D2hm3dq',...
%     'RLSplitAlphaForget_parameters_D1mCherry',...
%     'RLSplitAlphaForget_parameters_D2mCherry'};
% folder='Data/';
% names = {'D1hm4ditrialhistory10_4_18',...
%     'd2hm4ditrialhistory10_4_18',...
%     'd2hm3dqtrialhistory10_4_18',...
%     'D1mcherrytrialhistory10_4_18',...
%     'D2mcherrytrialhistory10_4_18'};
% names = {'D1hm4ditrialhistory10_4_18',...
%     'd2hm4ditrialhistory10_4_18',...
%     'd2hm3dqtrialhistory10_4_18',...
%     'D1mcherrytrialhistories_4.7.20',...
%     'D2mcherrytrialhistory10_4_18'};


AllData=[];
k=0;
for cohort=1:length(names)
    load([folder,names{cohort},'.mat'])
    for i=1:length(parameters.data)
        k=k+1;
        choices = parameters.data{i};
        entry = parameters.rawentrydata{i};
        exclude = parameters.rawdata{i}==0;
        entry(exclude)=[];
        rew = choices==1;
        disc=parameters.DiscTTC(i);
        rec=parameters.RecallTTC(i);
        phases = [ones(1,disc) 2*ones(1,rec)];
        impulsive = (phases'==2).*(entry);
        animal= k+0*phases;
        Group(k)=cohort;
        init=0*phases; init(1)=1;
        AllData = [AllData;[animal' phases' choices rew init' impulsive]];
        nt(k)=length(rew);
        ttc(k,:)=[disc,rec];
        
    end
end

nt = [nt sum(nt)];
ns=k;
ng=max(Group);
Choice=AllData(:,3)';
Reward=AllData(:,4)';
Subject=AllData(:,1)';
Init=AllData(:,5)';
Phase=AllData(:,2)';
PhaseChange = [0 Phase(2:end)~=Phase(1:end-1)];
PhaseChange(Init==1)=0;
Entry = AllData(:,6)';

%%

%Correct(18) = 4;
RL_data = struct('n_s',ns,'n_g',ng,'n_t',nt,'Choice',Choice,...
    'Reward',Reward,'Subject',Subject,'Group',Group,'Init',Init,'Phase',...
    Phase,'PhaseChange',PhaseChange,'Entry',Entry,'ttc',ttc);
end