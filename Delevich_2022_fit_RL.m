% Anne Collins 2022
% Code for Delevich et al 2022
% RL modeling
% this necessitates having Stan and MatlabStan installed


% to re-reun the fitting in stan, start at the top.
% to run only analyses on the fit model, start at cell #6 (load the fit)

%% 1. Load and preformat the data
clear all

RL_data = FormatData;

%% 2. set up and run the stan model

tic
% stan model parameters - set up for debugging, increase for running
nchains=2;%4;
niter = 250;%2500;
nw=100;%500;

fitRL = stan('file','RLaabb.stan','data',RL_data,'iter',niter,...
    'chains',nchains,'refresh',10,'warmup',nw,'thin',1);
fitRL.verbose = true;
%fitRL.check();% show progress
fitRL.block();% block further instructions

elapsed=toc;
disp('finished')
%save fitRL2 fitRL
%%%
summary=print(fitRL)
%% 3. extract samples and diagnostics
[samples,diagnostics] = extractsamples('matlabstan',fitRL);

% diagnostic report
%calculate convergence diagnostics based on the posterior samples
rtable = mcmctable(samples);
%print a report about all mcmc diagnostics
interpretdiagnostics(diagnostics,rtable)

%% 4. extract variables of interest and save them

samples=[];
for it=1:nchains
    mua_alpha_b1(it,:)=fitRL.sim.samples(it).mua_alpha_b1;
    mub_alpha_b1(it,:)=fitRL.sim.samples(it).mua_alpha_b1;
    mua_beta_b1(it,:)=fitRL.sim.samples(it).mua_beta_b1;
    mub_beta_b1(it,:)=fitRL.sim.samples(it).mub_beta_b1;
    mua_alpha_b2(it,:)=fitRL.sim.samples(it).mua_alpha_b2;
    mub_alpha_b2(it,:)=fitRL.sim.samples(it).mua_alpha_b2;
    mua_beta_b2(it,:)=fitRL.sim.samples(it).mua_beta_b2;
    mub_beta_b2(it,:)=fitRL.sim.samples(it).mub_beta_b2;
    mua_alpha_a1(it,:)=fitRL.sim.samples(it).mua_alpha_a1;
    mub_alpha_a1(it,:)=fitRL.sim.samples(it).mub_alpha_a1;
    mua_beta_a1(it,:)=fitRL.sim.samples(it).mua_beta_a1;
    mub_beta_a1(it,:)=fitRL.sim.samples(it).mub_beta_a1;
    mua_alpha_a2(it,:)=fitRL.sim.samples(it).mua_alpha_a2;
    mub_alpha_a2(it,:)=fitRL.sim.samples(it).mub_alpha_a2;
    mua_beta_a2(it,:)=fitRL.sim.samples(it).mua_beta_a2;
    mub_beta_a2(it,:)=fitRL.sim.samples(it).mub_beta_a2;
    alpha_b1_gr(it,:,:)=fitRL.sim.samples(it).alpha_b1_gr;
    beta_b1_gr(it,:,:)=fitRL.sim.samples(it).beta_b1_gr;
    alpha_b2_gr(it,:,:)=fitRL.sim.samples(it).alpha_b2_gr;
    beta_b2_gr(it,:,:)=fitRL.sim.samples(it).beta_b2_gr;
    alpha_a1_gr(it,:,:)=fitRL.sim.samples(it).alpha_a1_gr;
    beta_a1_gr(it,:,:)=fitRL.sim.samples(it).beta_a1_gr;
    alpha_a2_gr(it,:,:)=fitRL.sim.samples(it).alpha_a2_gr;
    beta_a2_gr(it,:,:)=fitRL.sim.samples(it).beta_a2_gr;
    a1(it,:,:)=fitRL.sim.samples(it).a1_ind;
    a2(it,:,:)=fitRL.sim.samples(it).a2_ind;
    b1s(it,:,:)=fitRL.sim.samples(it).b1_ind;
    b2s(it,:,:)=fitRL.sim.samples(it).b2_ind;
    Q0(it,:,:)=fitRL.sim.samples(it).Q0;
end
samples.mua_alpha_b1=[mua_alpha_b1(1,:) mua_alpha_b1(2,:) mua_alpha_b1(3,:)];
samples.mub_alpha_b1=[mub_alpha_b1(1,:) mub_alpha_b1(2,:) mub_alpha_b1(3,:)];
samples.mua_beta_b1=[mua_beta_b1(1,:) mua_beta_b1(2,:) mua_beta_b1(3,:)];
samples.mub_beta_b1=[mub_beta_b1(1,:) mub_beta_b1(2,:) mub_beta_b1(3,:)];
samples.mua_alpha_b2=[mua_alpha_b2(1,:) mua_alpha_b2(2,:) mua_alpha_b2(3,:)];
samples.mub_alpha_b2=[mub_alpha_b2(1,:) mub_alpha_b2(2,:) mub_alpha_b2(3,:)];
samples.mua_beta_b2=[mua_beta_b2(1,:) mua_beta_b2(2,:) mua_beta_b2(3,:)];
samples.mub_beta_b2=[mub_beta_b2(1,:) mub_beta_b2(2,:) mub_beta_b2(3,:)];
samples.mua_alpha_a1=[mua_alpha_a1(1,:) mua_alpha_a1(2,:) mua_alpha_a1(3,:)];
samples.mub_alpha_a1=[mub_alpha_a1(1,:) mub_alpha_a1(2,:) mub_alpha_a1(3,:)];
samples.mua_beta_a1=[mua_beta_a1(1,:) mua_beta_a1(2,:) mua_beta_a1(3,:)];
samples.mub_beta_a1=[mub_beta_a1(1,:) mub_beta_a1(2,:) mub_beta_a1(3,:)];
samples.mua_alpha_a2=[mua_alpha_a2(1,:) mua_alpha_a2(2,:) mua_alpha_a2(3,:)];
samples.mub_alpha_a2=[mub_alpha_a2(1,:) mub_alpha_a2(2,:) mub_alpha_a2(3,:)];
samples.mua_beta_a2=[mua_beta_a2(1,:) mua_beta_a2(2,:) mua_beta_a2(3,:)];
samples.mub_beta_a2=[mub_beta_a2(1,:) mub_beta_a2(2,:) mub_beta_a2(3,:)];
samples.a1=[squeeze(a1(1,:,:));squeeze(a1(2,:,:));squeeze(a1(3,:,:))];
samples.a2=[squeeze(a2(1,:,:));squeeze(a2(2,:,:));squeeze(a2(3,:,:))];
samples.b1s=[squeeze(b1s(1,:,:));squeeze(b1s(2,:,:));squeeze(b1s(3,:,:))];
samples.b2s=[squeeze(b2s(1,:,:));squeeze(b2s(2,:,:));squeeze(b2s(3,:,:))];
samples.alpha_b1_gr =[squeeze(alpha_b1_gr(1,:,:));squeeze(alpha_b1_gr(2,:,:));squeeze(alpha_b1_gr(3,:,:))];
samples.beta_b1_gr =[squeeze(beta_b1_gr(1,:,:));squeeze(beta_b1_gr(2,:,:));squeeze(beta_b1_gr(3,:,:))];
samples.mean_b1_gr=samples.alpha_b1_gr./samples.beta_b1_gr;
samples.alpha_b2_gr =[squeeze(alpha_b2_gr(1,:,:));squeeze(alpha_b2_gr(2,:,:));squeeze(alpha_b2_gr(3,:,:))];
samples.beta_b2_gr =[squeeze(beta_b2_gr(1,:,:));squeeze(beta_b2_gr(2,:,:));squeeze(beta_b2_gr(3,:,:))];
samples.mean_b2_gr=samples.alpha_b2_gr./samples.beta_b2_gr;
samples.alpha_a1_gr =[squeeze(alpha_a1_gr(1,:,:));squeeze(alpha_a1_gr(2,:,:));squeeze(alpha_a1_gr(3,:,:))];
samples.beta_a1_gr =[squeeze(beta_a1_gr(1,:,:));squeeze(beta_a1_gr(2,:,:));squeeze(beta_a1_gr(3,:,:))];
samples.mean_a1_gr=samples.alpha_a1_gr./samples.beta_a1_gr;
samples.alpha_a2_gr =[squeeze(alpha_a2_gr(1,:,:));squeeze(alpha_a2_gr(2,:,:));squeeze(alpha_a2_gr(3,:,:))];
samples.beta_a2_gr =[squeeze(beta_a2_gr(1,:,:));squeeze(beta_a2_gr(2,:,:));squeeze(beta_a2_gr(3,:,:))];
samples.mean_a2_gr=samples.alpha_a2_gr./samples.beta_a2_gr;
samples.Q0=[squeeze(Q0(1,:,:));squeeze(Q0(2,:,:));squeeze(Q0(3,:,:))];

%% 5. compute WAIC
% names = {'d1hm4di','d2hm4di',...
%     'd2hm3dq','mcherry'};
names = {'d1hm4di','d1mcherry','d2hm3dq','d2hm4di',...
   'd2mcherry','panhm4di'};
Q0=[];
Q0 = samples.Q0;
pars(:,:,1) = samples.b1s;
pars(:,:,2) = samples.b2s;
pars(:,:,3) = samples.a1;
pars(:,:,4) = samples.a2;
cd LLH_computation/
for k=1:size(Q0,1)
    [p,e,Q]=computeLLH_RL2ab(RL_data,Q0(k,:),squeeze(pars(k,:,:)));
    ps(k,:)=p;
    es(k,:,:)=e;
    Qs(k,:,:)=Q;
end
cd ..

lpd = sum(log(mean(ps)));
efpars=sum(var(log(ps)));
WAIC = -2*(lpd-efpars);

es=squeeze(mean(es));
Qs=squeeze(mean(Qs));
figure;hist(ps(:),100);

save('FitRLaabb2021', 'samples','RL_data', ...
    'WAIC','lpd','efpars','elapsed','summary','names','es','Qs','Group')
%% 6. load model

load('FitRLaabb2021', 'samples','RL_data', ...
    'WAIC','lpd','efpars','elapsed','summary','names','es','Qs')


%% 7. Q0

Q0=samples.Q0;
mean(Q0)
figure;
[N,X]=hist(Q0,30);
plot(X,N)
legend

%% 8. format group level parameters

pnames={'alpha1','alpha2','alphadiff','beta1','beta2','betadiff'};
grpars(:,:,1)=samples.mean_a1_gr;
grpars(:,:,2)=samples.mean_a2_gr;
grpars(:,:,3)=grpars(:,:,1)-grpars(:,:,2);
grpars(:,:,4)=samples.mean_b1_gr;
grpars(:,:,5)=samples.mean_b2_gr;
grpars(:,:,6)=grpars(:,:,4)-grpars(:,:,5);

figure
for p=1:size(grpars,3)
    subplot(2,3,p)
    [N,X]=hist(grpars(:,:,p),30);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',14)
    q=quantile(grpars(:,1,p)-grpars(:,2,p),[.025,.975]);
    [p q(1).*q(2)>0]
    title(pnames{p})
end
    legend(names)
% 
% diff = grpars(:,3,p)-grpars(:,4,p);
% figure;
% hist(diff,30)

%% 9. plot group parameters
clear ps
ps(1,:)=mean(samples.a1);
ps(2,:)=mean(samples.a2);
ps(3,:)=mean(samples.a1)-mean(samples.a2);
ps(4,:)=mean(samples.b1s);
ps(5,:)=mean(samples.b2s);
ps(6,:)=mean(samples.b2s)-mean(samples.b1s);
ps = ps';

figure
for p=1:6
    subplot(3,2,p)
    hold on
    for g=1:6
        T=find(RL_data.Group==g);
        errorbar(g,mean(ps(T,p)),std(ps(T,p))/sqrt(length(T)),'o','linewidth',2)
    end
    set(gca,'fontsize',14)
    title(pnames{p})
end
legend(names)
%% 10. format individual parameters
allps(:,:,1)=samples.a1;
allps(:,:,2)=samples.a2;
allps(:,:,3)=allps(:,:,1)-allps(:,:,2);
allps(:,:,4)=samples.b1s;
allps(:,:,5)=samples.b2s;
allps(:,:,6)=allps(:,:,4)-allps(:,:,5);
np=size(allps,3);
figure
for p=1:size(allps,3)
    for g=1:size(grpars,2)
        T=find(RL_data.Group==g);
        subplot(6,np,p+np*(g-1))
        hold on
        %plot(sort(allps(:,T,p)));
        [N,X]=hist([allps(:,T,p) grpars(:,g,p)],30);
        plot(X,N(:,1:end-1)./repmat(sum(N(:,1:end-1)),[size(N,1),1]),'linewidth',1);
        plot(X,N(:,end)/sum(N(:,end)),'k','linewidth',2);
        title(names{g})
        ylabel(pnames{p})
    end
end

%% 11. Plot group effects
for g=1:6
    T=find(RL_data.Group==g);
    meanps(:,g,:)=mean(allps(:,T,:),2);
end
figure
for p=1:size(meanps,3)
    subplot(3,size(meanps,3),p)
    [N,X]=hist(meanps(:,:,p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',16)
    disp(pnames{p})
    disp(['group ',names{1},' vs. ',names{2},': p, [CI95] = '])
    findpvalue(meanps(:,1,p)-meanps(:,2,p))
    disp(['group ',names{6},' vs. ',names{2},': p, [CI95] = '])
    findpvalue(meanps(:,6,p)-meanps(:,2,p))
    if p==size(meanps,3)
    legend(names)
    end
    title(pnames{p})
    
    subplot(3,size(meanps,3),size(meanps,3)+p)
    [N,X]=hist(meanps(:,[1 6],p)-meanps(:,[2 2],p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',14)
    title([pnames{p}, ' vs. d1mch'])
    legend(names([1 6]))
    if size(grpars,2)==6
        
    subplot(3,size(meanps,3),2*size(meanps,3)+p)
    [N,X]=hist(meanps(:,3:4,p)-meanps(:,[5 5],p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',14)
    title([pnames{p}, ' vs. d2mch'])
    legend(names([3:4]))
    disp(['group ',names{3},' vs. ',names{5},': p, [CI95] = '])
    [findpvalue(meanps(:,3,p)-meanps(:,5,p))]
    disp(['group ',names{4},' vs. ',names{5},': p, [CI95] = '])
    [findpvalue(meanps(:,4,p)-meanps(:,5,p))]
    disp(['group ',names{1},' vs. ',names{3},': p, [CI95] = '])
    [findpvalue(meanps(:,1,p)-meanps(:,3,p))]
    end
end

disp('parameter values group*parameter')
squeeze(mean(meanps))



%% 12. same - customize
for g=1:6
    T=find(RL_data.Group==g);
    meanps(:,g,:)=mean(allps(:,T,:),2);
end
figure
for p=1:size(meanps,3)
    subplot(3,size(meanps,3),p)
    [N,X]=hist(meanps(:,:,p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',16)
    disp(pnames{p})
    disp(['group ',names{1},' vs. ',names{2},': p, [CI95] = '])
    findpvalue(meanps(:,1,p)-meanps(:,2,p))
    disp(['group ',names{6},' vs. ',names{2},': p, [CI95] = '])
    findpvalue(meanps(:,6,p)-meanps(:,2,p))
    if p==size(meanps,3)
    legend(names)
    end
    title(pnames{p})
    
    subplot(3,size(meanps,3),size(meanps,3)+p)
    [N,X]=hist(meanps(:,[1 2 6],p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',14)
    title([pnames{p}, ' vs. d1mch'])
    legend(names([1 2 6]))
    if size(grpars,2)==6
        
    subplot(3,size(meanps,3),2*size(meanps,3)+p)
    [N,X]=hist(meanps(:,3:5,p),30);
    N = N./repmat(sum(N),[size(N,1),1]);
    plot(X,N,'-','linewidth',2);
    set(gca,'fontsize',14)
    title([pnames{p}, ' vs. d2mch'])
    legend(names([3:5]))
    disp(['group ',names{3},' vs. ',names{5},': p, [CI95] = '])
    [findpvalue(meanps(:,3,p)-meanps(:,5,p))]
    disp(['group ',names{4},' vs. ',names{5},': p, [CI95] = '])
    [findpvalue(meanps(:,4,p)-meanps(:,5,p))]
    disp(['group ',names{1},' vs. ',names{3},': p, [CI95] = '])
    [findpvalue(meanps(:,1,p)-meanps(:,3,p))]
    disp(['group ',names{1},' vs. ',names{6},': p, [CI95] = '])
    [findpvalue(meanps(:,1,p)-meanps(:,6,p))]
    end
end

squeeze(mean(meanps))

%% 13. validation - simulate

%RL_data=FormatData;
ttc=RL_data.ttc;

simttc=simulate_RL2ab(RL_data,mean(Q0),ps);

%% 14. plot validation
figure
for p=1:2
    %subplot(1,2,p)
    hold on
    plot(ttc(:,p),simttc(:,p),'o','markersize',6)
    lsline
    plot([min([ttc(:,p);simttc(:,1)]) max([ttc(:,p);simttc(:,1)])],...
        [min([ttc(:,p);simttc(:,1)]) max([ttc(:,p);simttc(:,1)])],'k')
    set(gca,'fontsize',14)
    xlabel('real ttc')
    ylabel(['ttc phase ',num2str(p)])
    title(['Phase ',num2str(p)])
end

%% 15. group validation figure
figure
for p=1:2
   for g=1:6
        T=find(RL_data.Group==g);
        subplot(2,2,p)
        plot(g+.1*randn(1,length(T)),ttc(T,p),'ok')
        hold on
        errorbar(g,mean(ttc(T,p)),std(ttc(T,p))/sqrt(length(T)),'+k','linewidth',2)
        set(gca,'fontsize',14, 'xtick',1:5,'xticklabel',names)
        ylabel('ttc')
        title('data')
        subplot(2,2,2+p)
        plot(g+.1*randn(1,length(T)),simttc(T,p),'ok')
        hold on
        errorbar(g,mean(simttc(T,p)),std(simttc(T,p))/sqrt(length(T)),'+k','linewidth',2)
        set(gca,'fontsize',14, 'xtick',1:6,'xticklabel',names)
        ylabel(['ttc phase ',num2str(p)])
        title('sim')
    end
end

%% 16. supplementary figure - group validation simulation
figure
for p=1:2
   for g=1:6
        T=find(RL_data.Group==g);
        subplot(1,2,p)
        hold on
        noise = g+.05*randn(length(T),1);
        plot(noise-.1,ttc(T,p),'or')
        plot(noise+.1,simttc(T,p),'*b')
        plot([noise-.1 noise+.1]',[ttc(T,p) simttc(T,p)]','-','color',[.5 .5 .5])
        errorbar(g-.1,mean(ttc(T,p)),std(ttc(T,p))/sqrt(length(T)),'ok','linewidth',2,'markeredgecolor','r')
        errorbar(g+.1,mean(simttc(T,p)),std(simttc(T,p))/sqrt(length(T)),'*k','linewidth',2,'markeredgecolor','b')
        plot([g-.1 g+.1],[mean(ttc(T,p)) mean(simttc(T,p))],'k','linewidth',2)
        set(gca,'fontsize',14, 'xtick',1:6,'xticklabel',names)
        ylabel(['ttc phase ',num2str(p)])
        title('sim')
    end
end


[rho,pval] = corr(ttc,simttc,'type','Spearman')