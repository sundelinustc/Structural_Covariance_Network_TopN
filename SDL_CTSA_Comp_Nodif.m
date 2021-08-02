function SDL_CTSA_Comp_Nodif(SDL)


% compare CT or SA between two groups 

%% load cleaned data
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'T');
% T    --- table containing values and sbj info
fprintf('Loaded: ComBat data loaded from <-%s\n\n\n',fn);

% %% lme model
% tbl = T(:,{'Group','Age','Gender','Site'}); 
% tbl.avg  = mean(T{:,2:149},2); % mean CT or SA for each participant
% k = 0;
% plist = []; % to contain the p value of effect of Group per area 
% R = []; % matrix containing residuals after lme
% for j = 2:149 % per column of data values  
%     k = k + 1;
%     tbl.data = T{:,j};
%     lme  = fitlme(tbl,'data ~ Group + Age + Age^2 + Gender + avg + (1|Site)'); % to test the Group effect
%     plist(k) = lme.Coefficients{2,6};
%     
%     lme1 = fitlme(tbl,'data ~         Age + Age^2 + Gender + avg + (1|Site)'); % to gain the residules
%     R(:,k) = residuals(lme1);
%     
%     fprintf('LME for Area No. %03d\n',k); 
% end
% fprintf('\n\nResults: Area No. showing significant (uncorrected) between-group difference in %s:\n',SDL.data_type{1});
% pl = find(plist<0.05)

%% T test
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT');
for j=1:148 % per area
    [~,p,~,stats] = ttest2(T{idx1,j+1},T{idx2,j+1}); % T columns of data are from 2 to 149
    d(j)          = computeCohen_d(T{idx1,j+1},T{idx2,j+1}); % Cohen's d, independent
    plist(j) = p;
    tlist(j) = stats.tstat;
end
find(plist<=0.05)

% FDR correction
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(plist,0.05,'pdep','yes');
fprintf('\n\n\n\n');


%% The Rank of Areas showing CT/SA Reduction
% [~,kk] = sort(p,'ascend'); % kk is the list of rank according to p-value (t test)
% [~,kk] = sort(d,'ascend'); % kk is the list of rank according to d-value (Cohen's d),from decreased CT/SA to increased CT/SA
[~,kk] = sort(abs(d),'ascend'); % kk is the list of rank according to absolute d-value (Cohen's d),from not much between-group differences in CT/SA to larger differences in CT/SA
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure('DefaultAxesFontSize',18)
barh(d(kk));
CC = strrep(T.Properties.VariableNames(kk(1:20:148)+1),'_', ' ');
yticklabels(CC);
xlabel('Effect Size');
title([SDL.data_type{1},' Changes in PTSD vs CONT']);
fdir = fullfile(SDL.out,SDL.data_type{1});
savefig(fullfile(fdir,['Figure Rank_Areas_Change_in_PTSD_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
fprintf('Completed: Ranked areas showing brain changes for PTSD vs CONT\n');

%% rank list of effect size
fdir = fullfile(SDL.path,'Original');
load(fullfile(fdir,'RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations
tbl = cell2table(D,'VariableNames',{'Area','X','Y','Z'});

fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Data_ComBat_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
T1 = load(fn); % load COmBat data for calculation of mean and SD


tbl.No(:)     = 1:148;
tbl.CohenD(:) = d;
tbl.tstat(:)  = tlist;
tbl.p(:)      = plist;
for i = 1:148 % per area
    tbl.PTSD(i)    = nanmean(T1.T{idx1,i+1}); % area No.in T is in fact from 2 to 149
    tbl.PTSDstd(i) = nanstd(T1.T{idx1,i+1});
    tbl.CONT(i)    = nanmean(T1.T{idx2,i+1});
    tbl.CONTstd(i) = nanstd(T1.T{idx2,i+1});
end

[~,I] = sort(kk); tbl.rank(:) = I;
tbl = sortrows(tbl,'rank'); % sort the table according to the rank

%% save results
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
save(fn,'tbl','d','kk');
fprintf('Saved: Effect Size (Cohen''s d) of %s reduction\n',SDL.data_type{1});

%% End
end