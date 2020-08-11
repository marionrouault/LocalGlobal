% Rouault & Fleming (2020) Formation of global self-beliefs in the human brain.
% Behavioral analyses and figures for the metacognition experiment (outside of the scanner).


function [] = BehaviorMetacogExp()

% ------------------------------------------------------------
%      Load data to reproduce plots and statistics
% ------------------------------------------------------------

% These are provided but can be re-generated using the HMeta-d toolbox
% (Fleming, 2017 Neurosci. Consciousness):
load fitHME
load fitHMD

load MetacogExp
load MainExp

% Accuracy, reaction times and confidence data:
trial_excluded = MetacogExp.trial_excluded ;
basic_perf     = MetacogExp.basic_perf ;
basic_conf     = MetacogExp.basic_conf ;
confid         = MetacogExp.confid     ;
bins_conf      = MetacogExp.bins_conf  ;

% Output from ordinal regression models predicting confidence:
conf_pred = MetacogExp.conf_pred   ;

% Number of participants:
nS = size(basic_perf,1);



% ---------------------------------------------------
%                To make the figures
% ---------------------------------------------------

% colors for graphs:
lightblue      = [.12 .42 .82];

bluecol        = [15 153 204]/255;
purplecol      = [153 15 204]/255;

colormodel1 = [.5 .2 .8];
colormodel2 = [.5 .8 .2];
colormodel3 = [25 235 187]/255;
colormodel4 = [204 15 153]/255;
colormodel5 = [.8 .5 .2];

% adding individual data points:
sss = .27 ;
ttt = .17 ;
frame = .1 ;



% ---------------------------------------------------
%           Make the statistical analyses
% ---------------------------------------------------

disp('=========== GROUP ANALYSIS ===========')
disp(['N = ',num2str(nS),' participants'])

disp(['Across subjects mean proportion of excluded trials: ', ...
    num2str(mean(trial_excluded))])

% Print basic performance and reaction times:
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

disp(['Performance on easy against difficult trials: tval = ',...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pv)])

[~,pva,~,stats] = ttest(basic_perf(:,3),basic_perf(:,4)) ;

disp(['RTs on easy against difficult trials: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pva)])


[~,pval,~,stats] = ttest(confid(:,1),confid(:,2)) ;
disp(['Confidence easy-corrects against difficult-corrects: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pval)])

[~,pvalu,~,stats] = ttest(confid(:,3),confid(:,4)) ;
disp(['Confidence easy-incorrects against diff-incorrects: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalu)])

[~,pvalue,~,stats] = ttest(basic_conf(:,2),basic_conf(:,3)) ;
disp(['Confidence incorrects against corrects: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalue)])

[~,pvalue,~,stats] = ttest(basic_conf(:,4),basic_conf(:,5)) ;
disp(['Confidence easy against diff. trials: tval = ', ...
    num2str(round(stats.tstat*100)/100),' pval = ',num2str(pvalue)])


% 2 x 2 ANOVA on confidence according to:
% Factor A [CORRECT,INCORRECT]
% Factor B [DIFF,EASY]
disp('ANOVA on confidence, factors Accuracy and Difficulty: ')
[~,F,~,p,~,~,~,~,~] = repanova(confid,[2 2]);




% ---------------------------------------------------
%                   Make the figures
% ---------------------------------------------------

% Objective performance and RTs for easy and difficult conditions:
figure(5)
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
title('Metacognition task','fontsize',25)
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


% Confidence according to accuracy:
figure(6)
hold on;
bar([mean(basic_conf(:,2)) 0],'FaceColor',[0 180 0]/255,'EdgeColor',[0 180 0]/255,'LineWidth',5)
bar([0 mean(basic_conf(:,3))],'FaceColor',[180 0 0]/255,'EdgeColor',[180 0 0]/255,'LineWidth',5)
for o = 1:nS
    hold on;
    plot([1+sss,1+ttt],[basic_conf(o,2),basic_conf(o,2)],'k-','LineWidth',2) ;
    plot([2+sss,2+ttt],[basic_conf(o,3),basic_conf(o,3)],'k-','LineWidth',2) ;
end
errorbar(1:2,mean(basic_conf(:,2:3)),std(basic_conf(:,2:3))/sqrt(nS),'k.','LineWidth',5)
ylabel('Confidence (raw ratings)','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Corrects','Errors',''})
axis([0 3 1 6])
hold off


% Confidence according to accuracy and difficulty level:
figure(3)
hold on;
errorbar(1:2,[mean(confid(:,1)) mean(confid(:,2))], ...
    [std(confid(:,1)) std(confid(:,2))]/sqrt(nS),'Color',[0 .6 0],'LineWidth',3)
errorbar(1:2,[mean(confid(:,3)) mean(confid(:,4))], ...
    [std(confid(:,3)) std(confid(:,4))]/sqrt(nS),'Color',[.6 0 0],'LineWidth',3)
ylabel('Confidence','fontsize',25)
set(gca,'fontsize',25,'LineWidth',1.5,'XTickLabel',{'','Diff','Easy',''})
axis([0 3 1 6])
hold off


% Confidence distributions:

% marginalise for the plot:
meanbins = mean(bins_conf,3) ;
stdbins  = std(bins_conf,0,3)/sqrt(nS) ;

numgroups = size(meanbins, 1);
numbars = size(meanbins, 2);

figure(4)
subplot(3,5,[6 7])
hold on;
h = bar(meanbins);
set(h,'BarWidth',1);
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, meanbins(:,i), stdbins(:,i),'k.','LineWidth',3)
end
legend({'Error','Correct'},'Location','northwest','Orientation','vertical')
axis([0 7 0 50])
xlabel('Confidence ratings','fontsize',15)
set(gca,'fontsize',15,'LineWidth',1.5,'XTicklabel',{'','1','2','3','4','5','6',''})
hold off




% ---------------------------------------------------
%       Metacognitive efficiency: H-Metad fits
% ---------------------------------------------------

% fit hierarchical HMeta-d model (Fleming, 2017 NOC toolbox):
%%% How the files provided were obtained:
% fitHME = fit_meta_d_mcmc_group(nR_S1_E, nR_S2_E) ;
% fitHMD = fit_meta_d_mcmc_group(nR_S1_D, nR_S2_D) ;

disp(['mean H-Mratio easy trials: ',num2str(exp(fitHME.mu_logMratio)), ...
    '. Convergence Rhat: ',num2str(fitHME.mcmc.Rhat.mu_logMratio)])
disp(['mean H-Mratio difficult trials: ',num2str(exp(fitHMD.mu_logMratio)), ...
    '. Convergence Rhat: ',num2str(fitHMD.mcmc.Rhat.mu_logMratio)])

plotSamples(exp(fitHME.mcmc.samples.mu_logMratio))
plotSamples(exp(fitHMD.mcmc.samples.mu_logMratio))

disp('HDI between easy and difficult trials:')
calc_HDI(exp(fitHME.mcmc.samples.mu_logMratio(:))-exp(fitHMD.mcmc.samples.mu_logMratio(:)))


% Metacognitive efficiency fitted separately for easy and difficult
% conditions:

% plot of posterior samples:
samplesE = exp(fitHME.mcmc.samples.mu_logMratio);
samplesD = exp(fitHMD.mcmc.samples.mu_logMratio);

figure(7)
% easy has a lighter color - hacked from plotSamples.m
hold on;
hist(samplesE(:),50);
set(get(gca,'child'),'FaceColor',lightblue,'EdgeColor','k')
hist(samplesD(:),50);
title('Metacognitive efficiency','fontsize',20)
xlabel('meta-d` / d`','fontsize',20)
ylabel('Posterior distributions','fontsize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
box off
hold off



% ----------------------------------------------------------------
%   Regressions for generalising confidence to the fMRI data
% ----------------------------------------------------------------

% Fitting statistics
betas_model1 = [];
betas_model2 = [];
betas_model3 = [];
betas_model4 = [];
betas_model5 = [];

dev_model1 = [];
dev_model2 = [];
dev_model3 = [];
dev_model4 = [];
dev_model5 = [];

% Sample betas corresponding only to regressors of interest:
for suj = 1:nS
    betas_model1 = [betas_model1 ; (conf_pred{1,suj}.STATSmodel1.beta(end-3:end))'];
    betas_model2 = [betas_model2 ; (conf_pred{1,suj}.STATSmodel2.beta(end-2:end))'];
    betas_model3 = [betas_model3 ; (conf_pred{1,suj}.STATSmodel3.beta(end-2:end))'];
    betas_model4 = [betas_model4 ; (conf_pred{1,suj}.STATSmodel4.beta(end-2:end))'];
    betas_model5 = [betas_model5 ; (conf_pred{1,suj}.STATSmodel5.beta(end-1:end))'];
    
    dev_model1 = [dev_model1 ; (conf_pred{1,suj}.DEVmodel1)];
    dev_model2 = [dev_model2 ; (conf_pred{1,suj}.DEVmodel2)];
    dev_model3 = [dev_model3 ; (conf_pred{1,suj}.DEVmodel3)];
    dev_model4 = [dev_model4 ; (conf_pred{1,suj}.DEVmodel4)];
    dev_model5 = [dev_model5 ; (conf_pred{1,suj}.DEVmodel5)];
end

% NB. the deviance of the fit at the solution vector. The deviance is a
% generalization of the residual sum of squares.

% Confidence according to accuracy and difficulty level:
% [confCorrDiff confCorrEasy confIncorrDiff confIncorrEasy]
% plotted on top of it: conf prediction from regression models 1-5:

confid_pred_model1 = zeros(nS,4) ;
confid_pred_model2 = zeros(nS,4) ;
confid_pred_model3 = zeros(nS,4) ;
confid_pred_model4 = zeros(nS,4) ;
confid_pred_model5 = zeros(nS,4) ;

for suj = 1:nS
    confid_pred_model1(suj,:) = conf_pred{suj}.confidmodel1 ;
    confid_pred_model2(suj,:) = conf_pred{suj}.confidmodel2 ;
    confid_pred_model3(suj,:) = conf_pred{suj}.confidmodel3 ;
    confid_pred_model4(suj,:) = conf_pred{suj}.confidmodel4 ;
    confid_pred_model5(suj,:) = conf_pred{suj}.confidmodel5 ;
end


figure(4)
subplot(3,5,8)
hold on;
errorbar(1:2,[mean(confid(:,1)) mean(confid(:,2))], ...
    [std(confid(:,1)) std(confid(:,2))]/sqrt(nS),'Color',[.6 0 0],'LineWidth',3)
errorbar(1:2,[mean(confid(:,3)) mean(confid(:,4))], ...
    [std(confid(:,3)) std(confid(:,4))]/sqrt(nS),'Color',[0 0 .6],'LineWidth',3)
errorbar(1:2+frame,[mean(confid_pred_model1(:,1)) mean(confid_pred_model1(:,2))], ...
    [std(confid_pred_model1(:,1)) std(confid_pred_model1(:,2))]/sqrt(nS),'Color',colormodel1,'LineWidth',3)
errorbar(1:2+frame,[mean(confid_pred_model1(:,3)) mean(confid_pred_model1(:,4))], ...
    [std(confid_pred_model1(:,3)) std(confid_pred_model1(:,4))]/sqrt(nS),'Color',colormodel1,'LineWidth',3)
ylabel('Predicted confidence','fontsize',15)
title('Metacognition task','fontsize',15)
set(gca,'fontsize',15,'LineWidth',1.5,'XTickLabel',{'','Diff.','Easy',''})
axis([0 3 2 5])
hold off


% Plot regression coefficients for all confidence models:

figure(4)
subplot(3,5,3)
hold on;
bar(mean(betas_model1),'FaceColor',[1 1 1],'EdgeColor',colormodel1,'LineWidth',3)
errorbar(1:size(betas_model1,2),mean(betas_model1),std(betas_model1)/sqrt(nS),'k.','LineWidth',3)
xlabel('Predictors of confidence','fontsize',15)
title('Model 3','fontsize',15,'Color',colormodel1)
set(gca,'fontsize',8,'LineWidth',1.5,'XTickLabel',{'','ACC','RT','DIF','INT',''})
axis([0 5 -.8 .8])
hold off

subplot(3,5,5)
hold on;
bar(mean(betas_model2),'FaceColor',[1 1 1],'EdgeColor',colormodel2,'LineWidth',3)
errorbar(1:size(betas_model2,2),mean(betas_model2),std(betas_model2)/sqrt(nS),'k.','LineWidth',3)
set(gca,'fontsize',12,'LineWidth',1.5,'XTickLabel',{'','ACC','RT','DIF',''})
title('Model 5','fontsize',15,'Color',colormodel2)
axis([0 4 -.8 .8])
hold off

subplot(3,5,2)
hold on;
bar(mean(betas_model3),'FaceColor',[1 1 1],'EdgeColor',colormodel3,'LineWidth',3)
errorbar(1:size(betas_model3,2),mean(betas_model3),std(betas_model3)/sqrt(nS),'k.','LineWidth',3)
set(gca,'fontsize',12,'LineWidth',1.5,'XTickLabel',{'','ACC','DIF','INT',''})
title('Model 2','fontsize',15,'Color',colormodel3)
axis([0 4 -.8 .8])
hold off

subplot(3,5,4)
hold on;
bar(mean(betas_model4),'FaceColor',[1 1 1],'EdgeColor',colormodel4,'LineWidth',3)
errorbar(1:size(betas_model4,2),mean(betas_model4),std(betas_model4)/sqrt(nS),'k.','LineWidth',3)
set(gca,'fontsize',12,'LineWidth',1.5,'XTickLabel',{'','ACC','RT','INT',''})
title('Model 4','fontsize',15,'Color',colormodel4)
axis([0 4 -.8 .8])
hold off

subplot(3,5,1)
hold on;
bar(mean(betas_model5),'FaceColor',[1 1 1],'EdgeColor',colormodel5,'LineWidth',3)
errorbar(1:size(betas_model5,2),mean(betas_model5),std(betas_model5)/sqrt(nS),'k.','LineWidth',3)
ylabel('Regression weights (a.u.)','fontsize',15)
set(gca,'fontsize',15,'LineWidth',1.5,'XTickLabel',{'','ACC','RT',''})
title('Model 1','fontsize',15,'Color',colormodel5)
axis([0 3 -.8 .8])
hold off

end
