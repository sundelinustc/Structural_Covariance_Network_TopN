function SDL_SCA_All_PTSDsev(SDL)
% Calculate the between-group differences of SCA
% also corect for multiple comparisons using FDR


%% Loading information
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = T.PTSDsev==1; % high severity in PTSD
idx2 = T.PTSDsev==0; % low severity in PTSD

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
S0 = SDL_SCA(R,idx1,idx2,func,0);     % SC of original PTSD & CONT data
S1 = SDL_SCA(R,idx1,idx2,func,SDL.N);  % SC of N times' permutations of PTSD & CONT labels, for whole-brain betwen-group contrast
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
% S0.M1 = [];S0.M2 = []; S1.M1 = [];S1.M2 = []; % to save space cause Matlab has difficulty to save data larger than 2GB
save(fn,'S0','S1','kk','-v7.3');
fprintf('Completed: Structural covariance matrix saved in -> %s\n\n\n',fn);




%% Load data
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn);
fprintf('Completed: Structural covariance matrix loaded from <- %s\n\n\n',fn);

% between-group difference
S0.Y1Y2 = triu(S0.Y1-S0.Y2,1); % original values of group1 - group 2, only showing values above diagnose
for ii = 1:size(S1.Y1,3) % per permutation
    S1.Y1Y2(:,:,ii) = triu(S1.Y1(:,:,ii)-S1.Y2(:,:,ii),1); % permuted values of group1 - group 2, only showing values above diagnose
end

% 95% confidence interval (CI)
diffCI = zeros(size(S1.Y1,1),size(S1.Y1,2),2);% :,:,1 - CI lower limit, :,:,2 - CI upper limit
pval   = ones(size(S1.Y1,1),size(S1.Y1,2));
for ii = 1:size(S1.Y1,1)
    for jj = 1:size(S1.Y1,2)
        if ii < jj % only above the diagnose
            v0 = S0.Y1Y2(ii,jj); % original between-group difference
            v1 = reshape(S1.Y1Y2(ii,jj,:),1,size(S1.Y1,3)); % permuted between-group difference
            dCI = SDL_CI(v1',0.95);
            diffCI(ii,jj,1) = dCI(1); diffCI(ii,jj,2) = dCI(2);
            pval(ii,jj)   = SDL_p_permutation(v0,v1);
        end
    end
end

%% plot uncorrected p values;
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure();im=imagesc((1-pval).*(pval<0.05));im.AlphaData = 1;colormap('hot');colorbar; % highlight pvals < 0.05
savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Whole_PTSD_vs_CONT',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

%% corrected p values
pp = nonzeros(triu(pval,1)); % p values above the diagnose
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pp,0.05,'pdep','yes');
%!!!No significant outputs after FDR (whatever "pdep" or "dep") correction.

%% save results
fn = fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Whole_PTSD_vs_CONT',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'pval','adj_p','diffCI');
fprintf('Saved: Whole-brain analyses results saved in \n-->%s\n',fn);



%% End
end





