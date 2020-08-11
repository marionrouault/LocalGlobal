% Rouault & Fleming (2020) Formation of global self-beliefs in the human brain.
% Behavioral analyses and figures for the main fMRI experiment.


function [] = BehaviorMainExp()


% ---------------------------------------------------
%        Load data to reproduce plots
% ---------------------------------------------------

load MainExp MainExp

% Objective performance in learning block trials:
trial_excluded  = MainExp.trial_excluded  ;
meanAccLB       = MainExp.meanAccLB       ;
basic_perf      = MainExp.basic_perf      ;
acc_progression = MainExp.acc_progression ;

% Global self-performance estimates.
% Task choice as a function of difference in accuracy between tasks:
x_accDiff_E     = MainExp.x_accDiff_E     ;
y_accDiff_E     = MainExp.y_accDiff_E     ;

% How often the easiest of both tasks was eventually chosen:
freqchE         = MainExp.freqchE         ;

% How often the best-performed of both tasks was eventually chosen:
freqchBest      = MainExp.freqchBest      ;

% Predicted confidence on scanner data using regressions from metacog task:
scan_confid_pred_model1 = MainExp.scan_confid_pred_model1 ;

% Task choice split as top vs. bottom of blocks
y_accDiff_E2    = MainExp.y_accDiff_E2    ;

% Link predicted confidence to main experiment task choices:
conf_C_I = MainExp.conf_C_I ;
delta_conf_beta = MainExp.delta_conf_beta ;
delta_acc_beta = MainExp.delta_acc_beta ;
delta_rt_beta = MainExp.delta_rt_beta ;
deviances = MainExp.deviances ;

% Regressions:
Xacc = MainExp.Xacc ;
Xdiff = MainExp.Xdiff ;
choseTask2 = MainExp.choseTask2 ;% coded 0=choseTask1, 1=choseTask2

% Number of participants
nS = length(freqchE);

% Number of learning blocks
N_series = 32 ;



% ------------------------------------------------------------------------
%                       Tools used for all figures
% ------------------------------------------------------------------------

% colors for graphs:
color_nofbeasy = [0 194 92]/255;
color_nofbdif  = [255 153 21]/255;
bluecol        = [15 153 204]/255;
purplecol      = [153 15 204]/255;
colormodel1 = [.5 .2 .8];

bleu_incong = [51 204 204]/255 ;
jaune_cong = [241 196 15]/255 ;

% adding individual data points:
sss = .27 ;
frame = .1 ;



% ---------------------------------------------------
%        Make the statistical analyses
% ---------------------------------------------------


disp('=========== GROUP ANALYSIS ===========')

disp(['N = ',num2str(nS),' participants'])

disp(['Across subjects mean proportion of excluded trials: ', ...
    num2str(mean(trial_excluded))])

disp(['Across subjects mean overall performance: ', ...
    num2str(mean(meanAccLB))])


% Basic performance and reaction times:
disp(['Across subjects mean performance in easy trials: ', ...
    num2str(mean(basic_perf(:,2)))])
disp(['Across subjects mean performance in difficult trials: ', ...
    num2str(mean(basic_perf(:,1)))])
disp(['Across subjects mean RTs in easy trials: ', ...
    num2str(mean(basic_perf(:,4)))])
disp(['Across subjects mean RTs in difficult trials: ', ...
    num2str(mean(basic_perf(:,3)))])


% Test for a difference between easy and difficult conditions:
[~,pv,~,stats] = ttest(basic_perf(:,1),basic_perf(:,2)) ;

disp(['Performance easy against difficult: tval = ',...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pv)])

[~,pva,~,stats] = ttest(basic_perf(:,3),basic_perf(:,4)) ;

disp(['RTs easy against difficult: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pva)])


disp(['Range of choosing the best-performed of both tasks: Min = ', ...
    num2str(min(freqchBest)*100),'%     Max = ', ...
    num2str(max(freqchBest)*100),'%'])


% Global self-performance estimates:
[~,pval,~,stats] = ttest(y_accDiff_E2(:,1),y_accDiff_E2(:,3)) ;
disp(['Frequency chose easy task, top vs. bottom of blocks: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pval)])


% Link predicted confidence to main experiment task choices:
[~,pvalu,~,stats] = ttest(delta_conf_beta) ;
disp(['One sample t-test for beta delta confidence predicting task choice: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalu)])

% Examine whether link is stronger than for accuracy or RTs:
[~,pvalu,~,stats] = ttest(delta_acc_beta) ;
disp(['One sample t-test for beta delta accuracy predicting task choice: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalu)])

[~,pvalu,~,stats] = ttest(delta_rt_beta) ;
disp(['One sample t-test for beta delta RTs predicting task choice: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalu)])

% Compare these three models:
disp 'sum of deviances for confidence, accuracy, RT models:'
sum(deviances)

% Correlation between confidence on chosen task and confidence on chosen
% minus unchosen global SPEs:
[rho,pval] = corrcoef(conf_C_I(:,1),conf_C_I(:,1)-conf_C_I(:,2));



% ---------------------------------------------------
%                Make the figures
% ---------------------------------------------------


% Objective performance and RTs in easy and difficult conditions
figure;
subplot(2,2,1)
hold on;
bar(mean(basic_perf(:,1:2)),'FaceColor',[1 1 1],'EdgeColor',bluecol,'LineWidth',5)
errorbar(1:2,mean(basic_perf(:,1:2)),std(basic_perf(:,1:2))/sqrt(nS),'k.','LineWidth',5)
for o = 1:nS
    hold on;
    plot(1+sss,basic_perf(o,1),'ko','LineWidth',2) ;
    plot(2+sss,basic_perf(o,2),'ko','LineWidth',2) ;
end
ylabel('Objective performance','fontsize',25)
title('Main fMRI expt','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Diff.','Easy',''})
axis([0 3 .5 1])
hold off

subplot(2,2,2)
hold on;
bar(mean(basic_perf(:,3:4)),'FaceColor',[1 1 1],'EdgeColor',purplecol,'LineWidth',5)
errorbar(1:2,mean(basic_perf(:,3:4)),std(basic_perf(:,3:4))/sqrt(nS),'k.','LineWidth',5)
for o = 1:nS
    hold on;
    plot(1+sss,basic_perf(o,3),'ko','LineWidth',2) ;
    plot(2+sss,basic_perf(o,4),'ko','LineWidth',2) ;
end
ylabel('Reaction times (ms)','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Diff.','Easy',''})
axis([0 3 500 1200])
hold off


% General evolution across the whole experimental session
subplot(2,2,[3 4])
hold on;
errorbar(1:N_series,mean(acc_progression),std(acc_progression)/sqrt(nS),'kx','LineWidth',5);
plot(1:N_series,mean(acc_progression),'k.');lsline
plot(1:N_series,.5*ones(1,N_series),'c--','LineWidth',2)
ylabel('Mean performance per block','fontsize',25)
xlabel('Learning blocks','fontsize',25)
title('General progression','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5)
axis([0 N_series+1 .4 1])
hold off



% Global self-performance estimates as assayed in end-of-block task choices

% Proportion of choice of easy over difficult tasks across blocks:
figure;
hold on;
bar(1:2,[mean(freqchE) 0],'FaceColor',color_nofbeasy,'EdgeColor',[1 1 1])
bar(1:2,[0 mean(1-freqchE)],'FaceColor',color_nofbdif,'EdgeColor',[1 1 1])
errorbar(1:2,[mean(freqchE) 0],[std(freqchE)/sqrt(nS) 0],'k.','LineWidth',2)
errorbar(1:2,[0 mean(1-freqchE)],[0 std(1-freqchE)/sqrt(nS)],'k.','LineWidth',2)
for o = 1:nS
    hold on;
    plot(1+sss,freqchE(o),'ko','LineWidth',2) ;
    plot(2+sss,1-freqchE(o),'ko','LineWidth',2) ;
end
set(gca,'fontsize',22,'LineWidth',1.5,'XTickLabel',{'','Easy','Diff.',''})
title('Global self-performance estimates','fontsize',22)
ylabel('Task choice frequency','fontsize',22)
axis([0 length(1:2)+1 0 1])
hold off


figure;
hold on;
bar([mean(conf_C_I(:,1)) 0],'EdgeColor',[1 1 1],'FaceColor',jaune_cong,'LineWidth',5)
bar([0 mean(conf_C_I(:,2))],'EdgeColor',[1 1 1],'FaceColor',bleu_incong,'LineWidth',5)
errorbar(1:2,mean(conf_C_I),std(conf_C_I)/sqrt(nS),'k.','LineWidth',5)
ylabel('Predicted confidence','fontsize',25)
title('Global self-performance estimates','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Higher','Lower',''})
axis([0 3 3 4])
hold off


% Task choice split according to fluctuations in performance across blocks:

% filter subjects at ceiling on easy trials (no empty cell): sub 34 & 36
y_accDiff_E_plot = y_accDiff_E([1:33 35 37:41],:);
x_accDiff_E_plot = x_accDiff_E([1:33 35 37:41],:);

figure;
hold on;
errorbar(nanmean(x_accDiff_E_plot),nanmean(y_accDiff_E_plot),nanstd(y_accDiff_E_plot)/sqrt(nS-2),...
    'Color',color_nofbeasy,'LineStyle','-','Marker','o','MarkerEdgeColor',...
    color_nofbeasy,'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',3)
xlabel('Easy-Diff. task performance','fontsize',25)
ylabel('Easy task choice proportion','fontsize',25)
set(gca,'fontsize',22,'LineWidth',1.5)
axis([-.5 .5 0.2 0.8])
hold off

[h,p,~,stats]=ttest(y_accDiff_E_plot(:,1),y_accDiff_E_plot(:,3));



% ------------------------------------------------------------------------
%       Predicted confidence from ordinal regressions in Metacog task
% ------------------------------------------------------------------------

% NB: Numerotation of Model(i) different from the paper.
figure;
hold on;
errorbar(1:2+frame,[mean(scan_confid_pred_model1(:,1)) mean(scan_confid_pred_model1(:,2))], ...
    [std(scan_confid_pred_model1(:,1)) std(scan_confid_pred_model1(:,2))]/sqrt(nS), ...
    'Color',colormodel1,'LineWidth',3)
errorbar(1:2+frame,[mean(scan_confid_pred_model1(:,3)) mean(scan_confid_pred_model1(:,4))], ...
    [std(scan_confid_pred_model1(:,3)) std(scan_confid_pred_model1(:,4))]/sqrt(nS), ...
    'Color',colormodel1,'LineWidth',3)
ylabel('Predicted confidence','fontsize',25)
title('Main fMRI expt','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Diff.','Easy',''})
axis([0 3 2 5])
hold off



% ---------------------------------------------------
%                Perform regressions
% ---------------------------------------------------

% Logistic regression examining if difference in performance between
% tasks and difficulty level and their interaction have an influence on
% task choices (in FFX):

% pool all regressors in the GLM function:
Xdiff(Xdiff==85) = 1 ;
Xdiff(Xdiff==70) = -1 ;

Xint  = Xacc.*Xdiff ;

% z-score regressors for commensurability of regression coefficients:
Xacc_  = (Xacc-mean(Xacc))/std(Xacc) ;
Xdiff_ = (Xdiff-mean(Xdiff))/std(Xdiff) ;
Xint_  = (Xint-mean(Xint))/std(Xint) ;

% Gram-Schmidt orthogonalization (optional)
% [Q,R] = qr_gs([Xacc_ Xdiff_ Xint_]) ;
% Xacc_  = Q(:,1) ;
% Xdiff_ = Q(:,2) ;
% Xint_  = Q(:,3) ;

X = [Xacc_ Xdiff_ Xint_] ;

% dependent variable
Y = boolean(choseTask2); % NB. coded 0=choseTask1, 1=choseTask2.

[betas,~,stats] = glmfit(X,Y,'binomial','link','logit');

disp('Logistic regression predicting task choice:')
disp(['Baseline should be NS: beta = ',num2str(betas(1)),' pval = ',num2str(stats.p(1))])
disp(['Influence of accuracy: beta = ',num2str(betas(2)),' pval = ',num2str(stats.p(2))])
disp(['Influence of diff. level: beta = ',num2str(betas(3)),' pval = ',num2str(stats.p(3))])
disp(['Interaction: beta = ',num2str(betas(4)),' pval = ',num2str(stats.p(4))])

end
