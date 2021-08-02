function SDL_SCA_All4(SDL)
% Calculate the between-group differences of SCA
% also corect for multiple comparisons using FDR


%% Loading information
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
% idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); 
if strcmp(SDL.data_type{1}(4:end),'Age')
    idx11 = strcmp(T.Group,'PTSD') & (T.Age<nanmedian(T.Age));  % group 11, PTSD & Age < 33 (median)
    idx12 = strcmp(T.Group,'PTSD') & (T.Age>=nanmedian(T.Age)); % group 12, PTSD & Age >= 33 (median)
    idx21 = strcmp(T.Group,'CONT') & (T.Age<nanmedian(T.Age));  % group 21, CONT & Age < 33 (median)
    idx22 = strcmp(T.Group,'CONT') & (T.Age>=nanmedian(T.Age)); % group 22, CONT & Age >= 33 (median)
elseif strcmp(SDL.data_type{1}(4:end),'Gender')
    idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'M');  % group 11, PTSD & Male
    idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'F');  % group 12, PTSD & Female
    idx21 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'M');  % group 21, CONT & Male
    idx22 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'F');  % group 22, CONT & Female
elseif strcmp(SDL.data_type{1}(4:end),'Dep')
%     try
%         idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'1');  % group 11, PTSD & depression
%         idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'0');  % group 12, PTSD & no-depression
%         idx21 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'1');  % group 21, CONT & depression
%         idx22 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'0');  % group 22, CONT & no-depression
%     catch
        idx11 = strcmp(T.Group,'PTSD') & (T.Dep==1);  % group 11, PTSD & depression
        idx12 = strcmp(T.Group,'PTSD') & (T.Dep==0);  % group 12, PTSD & no-depression
        idx21 = strcmp(T.Group,'CONT') & (T.Dep==1);  % group 21, CONT & depression
        idx22 = strcmp(T.Group,'CONT') & (T.Dep==0);  % group 22, CONT & no-depression
%     end
else
end

fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
% T    --- table containing raw CT values and sbj info
% R    --- residules of CT after lme regression
% idx1 --- row index of PTSD in T & R
% idx2 --- row index of CONT in T & R
% d    --- Cohen's d for PTSD vs CONT
% kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction  



%% Structural Covariance Matrix
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,0);     % SC of original PTSD & CONT data
S1 = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,SDL.N);  % SC of N times' permutations of PTSD & CONT labels, for whole-brain betwen-group contrast
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
% S0.M1 = [];S0.M2 = []; S1.M1 = [];S1.M2 = []; % to save space cause Matlab has difficulty to save data larger than 2GB
save(fn,'S0','S1','kk','-v7.3');
fprintf('Completed: Structural covariance matrix saved in -> %s\n\n\n',fn);




%% Load data
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn);
fprintf('Completed: Structural covariance matrix loaded from <- %s\n\n\n',fn);

% between-group difference
S0.Y11v21 = triu(S0.Y11-S0.Y21,1); % original values of group1 - group 2 for Age<33, only showing values above diagnose
S0.Y12v22 = triu(S0.Y12-S0.Y22,1); % original values of group1 - group 2 for Age>=33, only showing values above diagnose
for ii = 1:size(S1.Y11,3) % per permutation
    S1.Y11v21(:,:,ii) = triu(S1.Y11(:,:,ii)-S1.Y21(:,:,ii),1); % permuted values of group1 - group 2 for Age<33, only showing values above diagnose
    S1.Y12v22(:,:,ii) = triu(S1.Y12(:,:,ii)-S1.Y22(:,:,ii),1); % permuted values of group1 - group 2 for Age>=33, only showing values above diagnose
end

%% for Y11v21 (group 1 - group 2 for Age <33)
% 95% confidence interval (CI)
diffCI = zeros(size(S1.Y11,1),size(S1.Y11,2),2);% :,:,1 - CI lower limit, :,:,2 - CI upper limit
pval   = ones(size(S1.Y11,1),size(S1.Y11,2));
for ii = 1:size(S1.Y11,1)
    for jj = 1:size(S1.Y11,2)
        if ii < jj % only above the diagnose
            v0 = S0.Y11v21(ii,jj); % original between-group difference
            v1 = reshape(S1.Y11v21(ii,jj,:),1,size(S1.Y11,3)); % permuted between-group difference
            dCI = SDL_CI(v1',0.95);
            diffCI(ii,jj,1) = dCI(1); diffCI(ii,jj,2) = dCI(2);
            pval(ii,jj)   = SDL_p_permutation(v0,v1);
        end
    end
end
% plot uncorrected p values
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure();im=imagesc((1-pval).*(pval<0.05));im.AlphaData = 1;colormap('hot');colorbar; % highlight pvals < 0.05
savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Whole_PTSD_vs_CONT_1',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
% corrected p values
pp = nonzeros(triu(pval,1)); % p values above the diagnose
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pp,0.05,'pdep','yes');
pval_Y11v21 = pval;
adj_p_Y11v21 = adj_p;
diffCI_Y11v21 = diffCI;


%% for Y12v22 (group 1 - group 2 for Age >=33)
% 95% confidence interval (CI)
diffCI = zeros(size(S1.Y12,1),size(S1.Y12,2),2);% :,:,1 - CI lower limit, :,:,2 - CI upper limit
pval   = ones(size(S1.Y12,1),size(S1.Y12,2));
for ii = 1:size(S1.Y12,1)
    for jj = 1:size(S1.Y12,2)
        if ii < jj % only above the diagnose
            v0 = S0.Y12v22(ii,jj); % original between-group difference
            v1 = reshape(S1.Y12v22(ii,jj,:),1,size(S1.Y12,3)); % permuted between-group difference
            dCI = SDL_CI(v1',0.95);
            diffCI(ii,jj,1) = dCI(1); diffCI(ii,jj,2) = dCI(2);
            pval(ii,jj)   = SDL_p_permutation(v0,v1);
        end
    end
end
% plot uncorrected p values
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure();im=imagesc((1-pval).*(pval<0.05));im.AlphaData = 1;colormap('hot');colorbar; % highlight pvals < 0.05
savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Whole_PTSD_vs_CONT_2',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
% corrected p values
pp = nonzeros(triu(pval,1)); % p values above the diagnose
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pp,0.05,'pdep','yes');
pval_Y12v22 = pval;
adj_p_Y12v22 = adj_p;
diffCI_Y12v22 = diffCI;


%% save results
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Whole_PTSD_vs_CONT',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'pval_Y11v21','adj_p_Y11v21','diffCI_Y11v21','pval_Y12v22','adj_p_Y12v22','diffCI_Y12v22');
fprintf('Saved: Whole-brain analyses results saved in \n-->%s\n',fn);



%% End
end





