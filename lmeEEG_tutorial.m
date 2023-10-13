%% lmeEEG: Tutorial
% 
% 
% Download input files:
% sEEG_all: Event-related EEG dataset simulated using the MATLAB-based toolbox SEREEGA (https://github.com/lrkrol/SEREEGA).  Simulated EEG data included a P1-N2-P3 complex with different intercepts for subjects (N = 30) and items (N = 10). Moreover, the P3 was differently modulated according to two experimental conditions (i.e., a two-level experimental factor: A vs. B). 
% chanlocs: channel location variable (EEGLAB format)
% sEEG_table: Event table specifying Subject (ID), Item, and experimental condition
websave('sEEG_all.mat','https://osf.io/download/xkj95/'); load('sEEG_all.mat'); 
websave('sEEG_table.mat', 'https://osf.io/download/8bkxz/'), load('sEEG_table.mat');
addpath(['..',filesep,'functions']);


%% STEP 1
% Conduct mixed models on each channel/timepoint combination.
ID = nominal(sEEG_table.ID); Item=nominal(sEEG_table.Item); CON = categorical(sEEG_table.Condition);
mEEG = nan(size(sEEG_all));
for ch = 1:size(sEEG_all,1)
    parfor tpoint = 1:size(sEEG_all,2)
        EEG = double(squeeze(sEEG_all(ch,tpoint,:)));
        EEG = table(EEG,CON,ID,Item);
        m = fitlme(EEG,'EEG~CON+(1|ID)+(1|Item)');
        mEEG(ch,tpoint,:) = fitted(m,'Conditional',0)+residuals(m); % Extract marginal EEG
    end
end

% Extract design matrix X
EEG = double(squeeze(sEEG_all(1,1,:)));
EEG = table(EEG,CON,ID,Item);
m = fitlme(EEG,'EEG~CON+(1|ID)+(1|Item)');
X = designMatrix(m);


%% STEP 2 
% Perform mass univariate linear regressions on “marginal” EEG data.
t_obs = nan(size(mEEG,1),size(mEEG,2),size(X,2));
betas = nan(size(mEEG,1),size(mEEG,2),size(X,2));
se = nan(size(mEEG,1),size(mEEG,2),size(X,2));
for ch = 1:size(mEEG,1)
    parfor tpoint = 1:size(mEEG,2)
        EEG = squeeze(mEEG(ch,tpoint,:));
        [t_obs(ch,tpoint,:), betas(ch,tpoint,:), se(ch,tpoint,:)]=lmeEEG_regress(EEG,X)
    end
end


%% STEP 3
% Perform permutation testing ...
nperms=2000; % number of permutations
% [rperms] = lmeEEG_permutations(ID,nperms); % nperms within-subjects permutations of X (for stimuli-within-condition designs)
[rperms] = lmeEEG_permutations2(nperms, ID, Item); % within subjects and items permutations of X (for fully-crossed designs)
t_perms = nan(nperms,size(mEEG,1),size(mEEG,2),size(X,2)); % Initialize t-map
for p =1:nperms
    XX = X(rperms(:,p),:);
    for ch = 1:size(mEEG,1)
        parfor tpoint = 1:size(mEEG,2)
            EEG = squeeze(mEEG(ch,tpoint,:));
            [t_perms(p,ch,tpoint,:)]=lmeEEG_regress(EEG,XX);
        end
    end
end

%% ... and apply TFCE.
% This part requires ept_TFCE toolbox (https://github.com/Mensen/ept_TFCE-matlab)
for i = 2:size(X,2) % i from 2 since permutation testing of the intercept is not possible
    if ndims(t_obs) == 3 % channels x timepoints x fixed effects
        Results.(matlab.lang.makeValidName(m.CoefficientNames{i})) = lmeEEG_TFCE(squeeze(t_obs(:,:,i)),squeeze(t_perms(:,:,:,i)),chanlocs,[0.66 2]);
    elseif ndims(t_obs) == 4 % channels x frequencies x timepoints x fixed effects
        Results.(matlab.lang.makeValidName(m.CoefficientNames{i})) = lmeEEG_TFCE(squeeze(t_obs(:,:,:,i)),squeeze(t_perms(:,:,:,:,i)),chanlocs,[0.66 2]);
    end
end


%% PLOT RESULTS
mT = Results.CON_ConditionB.Obs;
mT2 = mT;
mT(not(Results.CON_ConditionB.Mask))=0;

intercept = squeeze(betas(:,:,1));
beta = squeeze(betas(:,:,2));

% Raster plot
mT=mT(:,sEEG_times>=0);
tick_labels = struct2table(chanlocs).labels;
figure,
imagesc(mT)
xlim([0 101])
set(gca,'ytick',1:19,'FontSize',8,'FontName','Arial');
set(gca,'TickLength',[0 0]);
set(gca,'XTick',linspace(1,101,6),'XTickLabel',0:200:1000,'FontSize',8,'FontName','Arial');
cmap2=cmap; cmap2(129,:)=[.8 .8 .8];
set(gca,'clim',[-20 20],'colormap',cmap2)
yticklabels(tick_labels);
xlabel("Time (ms)","FontWeight","bold","FontSize",10, "FontName","Arial");
set(gca,'color','none')
hc=colorbar;
a = get(hc,'YTickLabel');
set(hc, 'colormap', cmap,'YTickLabel',a,'FontSize',8,'FontName','Arial');
ylabel(hc,'t-value','FontWeight','bold','FontSize',10,'FontName','Arial');

% TOPOPLOT
[~, chm] = ismember({'Fz' 'F3' 'F4' 'Cz'}, struct2table(chanlocs).labels);
figure,
eeglab nogui
cmap =colormap('jet');
topoplot(mean(mT2(:,38:44),2), chanlocs, 'style', 'map', 'gridscale', 300, ...
    'maplimits', [-20 20], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{chm,'o','k',4,1},'colormap',cmap)

%TRACEPLOT
figure,
y1 = mean(intercept(chm,:),1);
y2 = mean(beta(chm,:),1);
y3 = mean(intercept(chm,:)+beta(chm,:),1);
plot(sEEG_times,y1,"Color",[30 136 229]/255,"LineWidth",3);
hold on, plot(sEEG_times,y2,"Color",[216 27 96]/255,"LineWidth",3);
hold on, plot(sEEG_times,y3,"Color",[255 193 7]/255,"LineWidth",3);
yline(0), xline(0), xlim([-100 1000])
a = get(gca,"XTickLabel"); set(gca,"XTickLabel",a,"FontSize",8,"FontName","Arial");
a = get(gca,"YTickLabel"); set(gca,"YTickLabel",a,"FontSize",8,"FontName","Arial");
xlabel("Time (ms)","FontWeight","bold","FontSize",10, "FontName","Arial");
ylabel("Beta value (μV)","FontWeight","bold","FontSize",10, "FontName","Arial");
legend("β_0 (Intercept)", "β_1 (Condition B vs. Condition A)", "β_0 + β_1", Location="northeast");
