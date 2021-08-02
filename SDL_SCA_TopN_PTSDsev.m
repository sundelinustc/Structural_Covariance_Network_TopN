function SDL_SCA_TopN_PTSDsev(SDL)

% SCA within the top-N areas with CT or SA reduction


%% Loading information
% Residuals
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = T.PTSDsev==1; % high severity in PTSD
idx2 = T.PTSDsev==0; % low severity in PTSD

% Effect size
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'kk','d'); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
% T    --- table containing raw CT values and sbj info
% R    --- residules of CT after lme regression
% idx1 --- row index of PTSD in T & R
% idx2 --- row index of CONT in T & R
% d    --- Cohen's d for PTSD vs CONT
% kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction  

% Region's label & XYZ coordinates
fdir = fullfile(SDL.path,'Original');
load(fullfile(fdir,'RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations
tbl = cell2table(D,'VariableNames',{'Area','X','Y','Z'});


%% Structural Covariance Matrix -- Original
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA(R,idx1,idx2,func,0);     % SC of original PTSD & CONT data
% S0.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
% S0.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2


%% Generate multiple sets of random nodes with matching hemispherical distribution & distance
SC0 = []; SC1 = []; % mean SC values of top-N, original(0) & permutation(1) data
NR = SDL.N; % number of sets of N random nodes, default=5000
for i = 2:length(kk) % per top-N
    % original nodes
    SC0(i).node   = kk(1:i); % No. of the top-N areas
    
    % random nodes
    [jlist,md,dlist] = SDL_Rand(SC0(i).node,tbl,NR); % the No. of the top-N areas
    
    % matching Euclidian distance
    % it is weired that can't use ttest(dlist-md);
    k = 0;
    while (k < SDL.Nr) && (signrank(dlist-md) < 0.05) % loop<3000 times if the mean distance of random nodes is different from original nodes distance
        k = k + 1;
        if mean(dlist) < md % if mean distance of random nodes is smaller than md
            d1 = mean(dlist);
            while d1 < md
                [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
            end
            idx = find(dlist==min(dlist)); % replace the set with minimal distance with a new set with distance > md
            idx = idx(1);
            jlist(idx,:) = j1;
            dlist(idx) = d1;
        elseif mean(dlist) > md % if mean distance of random nodes is larger than md
            d1 = mean(dlist);
            while d1 > md
                [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
            end
            idx = find(dlist==max(dlist)); % replace the set with maximal distance with a new set with distance < md
            idx = idx(1);
            jlist(idx,:) = j1;
            dlist(idx) = d1;
        else
        end
    end
    
    SC1(i).node   = jlist; % No. of the top-N areas
    fprintf('Generated: randomly chosen %d sets of %3d nodes with matching hemisphere and distance\n',NR,i);
    
end
%save
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_PTSD_HighvsLow_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'S0','SC0','SC1','-v7.3');
fprintf('Saved: Top-N analyses midway results saved in \n-->%s\n',fn);

%% SC, CI & p values of top-N areas with CT/SA reduction
% loading
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_PTSD_HighvsLow_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'S0','SC0','SC1');
fprintf('Loaded: Top-N analyses midway results loaded from \n<--%s\n',fn);

% mean SC, CI and p values
CIType = 0.95; % 95% confidene interval (CI) 
for i = 2:size(SC0,2) % per top-N
    % original data
    node0 = SC0(i).node;
    SC0(i).Z1 = triu(S0.Y1(node0,node0),1); SC0(i).Z2 = triu(S0.Y2(node0,node0),1); % original data, values above diagonal
    SC0(i).Z1c = nonzeros(SC0(i).Z1); SC0(i).Z2c = nonzeros(SC0(i).Z2); % non-zero values, Num of connections X 1
    SC0(i).mean1 = mean(SC0(i).Z1c); % mean of all connections, 1X1, original data, group 1
    SC0(i).mean2 = mean(SC0(i).Z2c); % mean of all connections, 1X1, original data, group 2
    SC0(i).Z1 = []; SC0(i).Z2 = [];SC0(i).Z1c = []; SC0(i).Z2c = [];
    
    % permuted data
    for j = 1:size(SC1(i).node,1) % per permutation (totally 5000 times)
        node1 = SC1(i).node(j,:);
        SC1(i).Z1(:,:,j) = triu(S0.Y1(node1,node1),1); SC1(i).Z2(:,:,j) = triu(S0.Y2(node1,node1),1); % original data, values above diagonal
        SC1(i).Z1c(:,j)  = nonzeros(SC1(i).Z1(:,:,j)); SC1(i).Z2c(:,j) = nonzeros(SC1(i).Z2(:,:,j)); % non-zero values, Num of connections X Num of permutations
    end
    SC1(i).mean1 = mean(SC1(i).Z1c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 1
    SC1(i).mean2 = mean(SC1(i).Z2c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 2
    SC1(i).Z1 = []; SC1(i).Z2 = [];SC1(i).Z1c = []; SC1(i).Z2c = [];
    
    % confidence interval (CI)
    SC1(i).CI1(:)    = SDL_CI(SC1(i).mean1',CIType); % 95% CI
    SC1(i).CI2(:)    = SDL_CI(SC1(i).mean2',CIType); % 95% CI
    SC1(i).diffCI(:) = SDL_CI(SC1(i).mean1' - SC1(i).mean2',CIType); % 95% CI
    
    % p values
    SC1(i).p1     = SDL_p_permutation(SC0(i).mean1,SC1(i).mean1); % group1 vs permutation
    SC1(i).p2     = SDL_p_permutation(SC0(i).mean2,SC1(i).mean2); % group2 vs permutation
    SC1(i).p_diff = SDL_p_permutation(SC0(i).mean1-SC0(i).mean2,SC1(i).mean1-SC1(i).mean2); % group1-2 vs permutation
    
    fprintf('Calculated: SC, CI & p-val for Top-%3d Areas\n',i);
end

%save
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_HighvsLow_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'SC0','SC1','-v7.3');
fprintf('Saved: Top-N analyses results saved in \n-->%s\n',fn);


%% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT 
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_HighvsLow_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    X1(i-1) = SC0(i).mean1;
    X2(i-1) = SC0(i).mean2;
    p1(i-1) = SC1(i).p1;
    p2(i-1) = SC1(i).p2;
    
    X1R(i-1) = mean(SC1(i).mean1);
    X2R(i-1) = mean(SC1(i).mean2);
    
    X1RA(i-1,:) = SC1(i).mean1;
    X2RA(i-1,:) = SC1(i).mean2;

    CI1L(i-1) = SC1(i).CI1(1);
    CI1R(i-1) = SC1(i).CI1(2);
    
    CI2L(i-1) = SC1(i).CI2(1);
    CI2R(i-1) = SC1(i).CI2(2);
    
    diffL(i-1) = SC1(i).diffCI(1);
    diffR(i-1) = SC1(i).diffCI(2);
    p_diff(i-1) = SC1(i).p_diff;
end


%% FDR corrected p values per top-area
fprintf('\n\nFDR correction: PTSD High vs random\n')
[~, ~, ~, adj_p1] = fdr_bh(p1,0.05,'pdep','yes'); 
find(adj_p1<=0.05)
fprintf('\n\nFDR correction: PTSD Low vs random\n')
[~, ~, ~, adj_p2] = fdr_bh(p2,0.05,'pdep','yes');
find(adj_p2<=0.05)
fprintf('\n\nFDR correction: PTSD High-Low vs random\n')
[~, ~, ~, adj_p_diff] = fdr_bh(p_diff,0.05,'pdep','yes');
find(adj_p_diff<=0.05)


%% p value for area under curve
p1 = SDL_p_permutation(mean(X1),mean(X1RA,1))
p2 = SDL_p_permutation(mean(X2),mean(X2RA,1))
if mean(X1)-mean(X2) > 0
    p_diff = SDL_p_permutation(mean(X1)-mean(X2),mean(X1RA,1)-mean(X2RA,1))
else
    p_diff = SDL_p_permutation(mean(X2)-mean(X1),mean(X2RA,1)-mean(X1RA,1))
end


%% Plot curves & heatmap of SC values
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure; % original vs permutation
subplot(3,2,1);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI1L, fliplr(CI1R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X1,'r-',X,X1R,'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD High');

subplot(3,2,3);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI2L, fliplr(CI2R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X2,'r-',X,X2R,'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD Low');

subplot(3,2,5);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [diffL, fliplr(diffR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,(X1-X2),'r-',X,(X1R-X2R),'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD High-Low');

MN1 = []; MN2 = []; MNd = [];
% re-organize matrix based on kk
for i = 1:length(kk)
    for j = 1:length(kk)
        MN1(i,j) = S0.Y1(kk(i),kk(j));
        MN2(i,j) = S0.Y2(kk(i),kk(j));
    end
end
MNd = MN1 - MN2;
% figure; 
subplot(3,2,2);imagesc(MN1); colormap('hot'); colorbar; title('PTSD High'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
subplot(3,2,4);imagesc(MN2); colormap('hot'); colorbar; title('PTSD Low'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
subplot(3,2,6);imagesc(MNd); colormap('hot'); colorbar; title('PTSD High-Low'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([-0.1,0.1]);
savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_2Dmap_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));


%% Plot correlation
% Both pos and neg efect size
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure('DefaultAxesFontSize',18)
ax1 = subplot(1,2,1); % PTSD
a = triu(S0.Y1,1); a = a + a'; % structural covariance after removing Inf in the diagnose
scatter(ax1,d,mean(a),45,'ko');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; axis([-0.13,0.13,-0.025,0.025]);
title('PTSD High'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
[r,p] = corr(d',mean(a)');
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(0,0.015,txt,'FontSize',18);

ax2 = subplot(1,2,2); % CONT
a = triu(S0.Y2,1); a = a + a'; % structural covariance after removing Inf in the diagnose
scatter(ax2,d,mean(a),45,'ko');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  axis([-0.13,0.13,-0.025,0.025]);
title('PTSD Low'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
[r,p] = corr(d',mean(a)');
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(0,0.015,txt,'FontSize',18);

savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Corr_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

% Only neg efect size
figure('DefaultAxesFontSize',18)
ax1 = subplot(1,2,1); % PTSD
a = triu(S0.Y1,1); a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',mean(a)']; a = a(a(:,1)<0,:);
scatter(ax1,a(:,1),a(:,2),45,'ko');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; axis([-0.11,0.01,-0.025,0.025]);
title('PTSD High'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);

ax2 = subplot(1,2,2); % CONT
a = triu(S0.Y2,1); a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',mean(a)']; a = a(a(:,1)<0,:);
scatter(ax2,a(:,1),a(:,2),45,'ko');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  axis([-0.11,0.01,-0.025,0.025]);
title('PTSD Low'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);

savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Corr_NegEffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

%% End
end