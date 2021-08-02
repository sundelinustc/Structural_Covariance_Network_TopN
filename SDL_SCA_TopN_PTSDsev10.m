function SDL_SCA_TopN_PTSDsev10(SDL)

% SCA within the top-N areas with CT or SA reduction


%% Loading information
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
T = T(T.CAPStype==4,:); % only subjects with CAPS-IV are included
R = T{:,2:149};


a = prctile(T.CAPS,[20 40 60 80]) % divide into 5 groups with compariable sample size
a = [4 20 46 68] % based on SA data, CT data are similiar
[sum(T.CAPS<a(1)),...
    sum(T.CAPS>=a(1)&T.CAPS<a(2)),...
    sum(T.CAPS>=a(2)&T.CAPS<a(3)),...
    sum(T.CAPS>=a(3)&T.CAPS<a(4)),...
    sum(T.CAPS>=a(4))]

idx1 = T.CAPS< a(1)               & T.CAPStype==4;   % group 1
idx2 = T.CAPS>=a(1) & T.CAPS<a(2) & T.CAPStype==4;   % group 2
idx3 = T.CAPS>=a(2) & T.CAPS<a(3) & T.CAPStype==4;   % group 3
idx4 = T.CAPS>=a(3) & T.CAPS<a(4) & T.CAPStype==4;   % group 4
idx5 = T.CAPS>=a(4) & T.CAPS<=120 & T.CAPStype==4;   % group 5

% idx1 = T.CAPS<2;                 % group 1
% idx2 = T.CAPS>=2  & T.CAPS<18;   % group 2
% idx3 = T.CAPS>=18 & T.CAPS<44;   % group 3
% idx4 = T.CAPS>=44 & T.CAPS<67;   % group 4
% idx5 = T.CAPS>=67 & T.CAPS<=120; % group 5

fdir = fullfile(SDL.out,SDL.data_type{1}(1:2));
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
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
S0 = SDL_SCA10_sev(R,idx1,idx2,idx3,idx4,idx5,func,0);     % SC of original data

% %% Generate multiple sets of random nodes with matching hemispherical distribution & distance
% SC0 = []; SC1 = []; % mean SC values of top-N, original(0) & permutation(1) data
% NR = 5000; % number of sets of N random nodes
% for i = 2:length(kk) % per top-N
%     % original nodes
%     SC0(i).node   = kk(1:i); % No. of the top-N areas
%     
%     % random nodes
%     [jlist,md,dlist] = SDL_Rand(SC0(i).node,tbl,NR); % the No. of the top-N areas
%     
%     % matching Euclidian distance
%     % it is weired that can't use ttest(dlist-md);
%     k = 0;
%     while (k < 3000) && (signrank(dlist-md) < 0.05) % loop<3000 times if the mean distance of random nodes is different from original nodes distance
%         k = k + 1;
%         if mean(dlist) < md % if mean distance of random nodes is smaller than md
%             d1 = mean(dlist);
%             while d1 < md
%                 [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
%             end
%             idx = find(dlist==min(dlist)); % replace the set with minimal distance with a new set with distance > md
%             idx = idx(1);
%             jlist(idx,:) = j1;
%             dlist(idx) = d1;
%         elseif mean(dlist) > md % if mean distance of random nodes is larger than md
%             d1 = mean(dlist);
%             while d1 > md
%                 [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
%             end
%             idx = find(dlist==max(dlist)); % replace the set with maximal distance with a new set with distance < md
%             idx = idx(1);
%             jlist(idx,:) = j1;
%             dlist(idx) = d1;
%         else
%         end
%     end
%     
%     SC1(i).node   = jlist; % No. of the top-N areas
%     fprintf('Generated: randomly chosen %d sets of %3d nodes with matching hemisphere and distance\n',NR,i);
%     
% end
% %save
% fdir = fullfile(SDL.out,SDL.data_type{1}); 
% fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% save(fn,'S0','SC0','SC1','-v7.3');
% fprintf('Saved: Top-N analyses midway results saved in \n-->%s\n',fn);

%% SC, CI & p values of top-N areas with CT/SA reduction
% loading
fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); % load SC0 and SC1 from existed results, i.e. CleanData, to avoid wasting time
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');
fprintf('Loaded: Top-N analyses midway results loaded from \n<--%s\n',fn);

% mean SC, CI and p values
CIType = 0.95; % 95% confidene interval (CI) 
for i = 2:size(SC0,2) % per top-N
    % original data, values above diagonal
    node0 = SC0(i).node;
    SC0(i).Z1 = triu(S0.Y1(node0,node0),1); 
    SC0(i).Z2 = triu(S0.Y2(node0,node0),1); 
    SC0(i).Z3 = triu(S0.Y3(node0,node0),1);
    SC0(i).Z4 = triu(S0.Y4(node0,node0),1);
    SC0(i).Z5 = triu(S0.Y5(node0,node0),1);
    
    
    % non-zero values, Num of connections X 1
    SC0(i).Z1c = nonzeros(SC0(i).Z1); 
    SC0(i).Z2c = nonzeros(SC0(i).Z2); 
    SC0(i).Z3c = nonzeros(SC0(i).Z3);
    SC0(i).Z4c = nonzeros(SC0(i).Z4);
    SC0(i).Z5c = nonzeros(SC0(i).Z5);

    
    % mean of all connections, 1X1, original data
    SC0(i).mean1raw = mean(SC0(i).Z1c); % mean of all connections, 1X1, original data, group 11
    SC0(i).mean2raw = mean(SC0(i).Z2c); % mean of all connections, 1X1, original data, group 12
    SC0(i).mean3raw = mean(SC0(i).Z3c);
    SC0(i).mean4raw = mean(SC0(i).Z4c);
    SC0(i).mean5raw = mean(SC0(i).Z5c);
    
    SC0(i).mean1pos = mean(SC0(i).Z1c(SC0(i).Z1c>0));
    SC0(i).mean2pos = mean(SC0(i).Z2c(SC0(i).Z2c>0));
    SC0(i).mean3pos = mean(SC0(i).Z3c(SC0(i).Z3c>0));
    SC0(i).mean4pos = mean(SC0(i).Z4c(SC0(i).Z4c>0));
    SC0(i).mean5pos = mean(SC0(i).Z5c(SC0(i).Z5c>0));
    
    SC0(i).mean1neg = mean(SC0(i).Z1c(SC0(i).Z1c<0));
    SC0(i).mean2neg = mean(SC0(i).Z2c(SC0(i).Z2c<0));
    SC0(i).mean3neg = mean(SC0(i).Z3c(SC0(i).Z3c<0));
    SC0(i).mean4neg = mean(SC0(i).Z4c(SC0(i).Z4c<0));
    SC0(i).mean5neg = mean(SC0(i).Z5c(SC0(i).Z5c<0));
    
    SC0(i).mean1abs = mean(abs(SC0(i).Z1c));
    SC0(i).mean2abs = mean(abs(SC0(i).Z2c));
    SC0(i).mean3abs = mean(abs(SC0(i).Z3c));
    SC0(i).mean4abs = mean(abs(SC0(i).Z4c));
    SC0(i).mean5abs = mean(abs(SC0(i).Z5c));
    
    
    SC0(i).Z1 = []; 
    SC0(i).Z2 = []; 
    SC0(i).Z3 = [];
    SC0(i).Z4 = [];
    SC0(i).Z5 = [];
    
    SC0(i).Z1c = []; 
    SC0(i).Z2c = [];
    SC0(i).Z3c = [];
    SC0(i).Z4c = [];
    SC0(i).Z5c = [];

    
    % permuted data
    for j = 1:size(SC1(i).node,1) % per permutation (totally 5000 times)
        node1 = SC1(i).node(j,:);
        SC1(i).Z1(:,:,j) = triu(S0.Y1(node1,node1),1); 
        SC1(i).Z2(:,:,j) = triu(S0.Y2(node1,node1),1); 
        SC1(i).Z3(:,:,j) = triu(S0.Y3(node1,node1),1);
        SC1(i).Z4(:,:,j) = triu(S0.Y4(node1,node1),1);
        SC1(i).Z5(:,:,j) = triu(S0.Y5(node1,node1),1);
        
        SC1(i).Z1c(:,j) = nonzeros(SC1(i).Z1(:,:,j)); 
        SC1(i).Z2c(:,j) = nonzeros(SC1(i).Z2(:,:,j)); 
        SC1(i).Z3c(:,j) = nonzeros(SC1(i).Z3(:,:,j));
        SC1(i).Z4c(:,j) = nonzeros(SC1(i).Z4(:,:,j));
        SC1(i).Z5c(:,j) = nonzeros(SC1(i).Z5(:,:,j));
        
        SC1(i).mean1pos(j) = mean(nonzeros(SC1(i).Z1c(:,j).*(SC1(i).Z1c(:,j)>0)));
        SC1(i).mean2pos(j) = mean(nonzeros(SC1(i).Z2c(:,j).*(SC1(i).Z2c(:,j)>0)));
        SC1(i).mean3pos(j) = mean(nonzeros(SC1(i).Z3c(:,j).*(SC1(i).Z3c(:,j)>0)));
        SC1(i).mean4pos(j) = mean(nonzeros(SC1(i).Z4c(:,j).*(SC1(i).Z4c(:,j)>0)));
        SC1(i).mean5pos(j) = mean(nonzeros(SC1(i).Z5c(:,j).*(SC1(i).Z5c(:,j)>0)));
        
        SC1(i).mean1neg(j) = mean(nonzeros(SC1(i).Z1c(:,j).*(SC1(i).Z1c(:,j)<0)));
        SC1(i).mean2neg(j) = mean(nonzeros(SC1(i).Z2c(:,j).*(SC1(i).Z2c(:,j)<0)));
        SC1(i).mean3neg(j) = mean(nonzeros(SC1(i).Z3c(:,j).*(SC1(i).Z3c(:,j)<0)));
        SC1(i).mean4neg(j) = mean(nonzeros(SC1(i).Z4c(:,j).*(SC1(i).Z4c(:,j)<0)));
        SC1(i).mean5neg(j) = mean(nonzeros(SC1(i).Z5c(:,j).*(SC1(i).Z5c(:,j)<0)));
    end
    SC1(i).mean1raw = mean(SC1(i).Z1c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 1
    SC1(i).mean2raw = mean(SC1(i).Z2c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 2
    SC1(i).mean3raw = mean(SC1(i).Z3c,1);
    SC1(i).mean4raw = mean(SC1(i).Z4c,1);
    SC1(i).mean5raw = mean(SC1(i).Z5c,1);
    
    SC1(i).mean1abs = mean(abs(SC1(i).Z1c),1); % mean of all absolute values of connections per permutation, 1 X Num of permutations, random data, group 1
    SC1(i).mean2abs = mean(abs(SC1(i).Z2c),1); % mean of all absolute values of connections per permutation, 1 X Num of permutations, random data, group 2
    SC1(i).mean3abs = mean(abs(SC1(i).Z3c),1);
    SC1(i).mean4abs = mean(abs(SC1(i).Z4c),1);
    SC1(i).mean5abs = mean(abs(SC1(i).Z5c),1);
    
    SC1(i).Z1 = []; 
    SC1(i).Z2 = [];
    SC1(i).Z3 = [];
    SC1(i).Z4 = [];
    SC1(i).Z5 = [];
    
    SC1(i).Z1c = []; 
    SC1(i).Z2c = [];
    SC1(i).Z3c = [];
    SC1(i).Z4c = [];
    SC1(i).Z5c = [];
    
%     % confidence interval (CI)
%     SC1(i).CI11(:)    = SDL_CI(SC1(i).mean11',CIType); % 95% CI
%     SC1(i).CI12(:)    = SDL_CI(SC1(i).mean12',CIType); % 95% CI
%     SC1(i).CI21(:)    = SDL_CI(SC1(i).mean21',CIType); % 95% CI
%     SC1(i).CI22(:)    = SDL_CI(SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_11v21(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean21',CIType); % 95% CI
%     SC1(i).diffCI_12v22(:) = SDL_CI(SC1(i).mean12' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_11v12(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean12',CIType); % 95% CI
%     SC1(i).diffCI_21v22(:) = SDL_CI(SC1(i).mean21' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_main1(:) = SDL_CI(SC1(i).mean11' + SC1(i).mean12' - SC1(i).mean21' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_main2(:) = SDL_CI(SC1(i).mean11' + SC1(i).mean21' - SC1(i).mean12' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_inter(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean12' - SC1(i).mean21' + SC1(i).mean22',CIType); % 95% CI
%     
%     % p values
%     SC1(i).p1     = SDL_p_permutation(SC0(i).mean1,SC1(i).mean1); % group1 vs permutation
%     SC1(i).p2     = SDL_p_permutation(SC0(i).mean2,SC1(i).mean2); % group2 vs permutation
%     SC1(i).p3     = SDL_p_permutation(SC0(i).mean3,SC1(i).mean3); % group3 vs permutation
%     SC1(i).p4     = SDL_p_permutation(SC0(i).mean4,SC1(i).mean4); % group4 vs permutation
%     SC1(i).p5     = SDL_p_permutation(SC0(i).mean5,SC1(i).mean5); % group5 vs permutation
% 
%     
%     SC1(i).p_diff_1v2 = SDL_p_permutation(SC0(i).mean1-SC0(i).mean2,SC1(i).mean1-SC1(i).mean2); % group 1-2 vs permutation
%     SC1(i).p_diff_1v3 = SDL_p_permutation(SC0(i).mean1-SC0(i).mean3,SC1(i).mean1-SC1(i).mean3); % group 1-3 vs permutation
%     SC1(i).p_diff_1v4 = SDL_p_permutation(SC0(i).mean1-SC0(i).mean4,SC1(i).mean1-SC1(i).mean4); % group 1-4 vs permutation
%     SC1(i).p_diff_1v5 = SDL_p_permutation(SC0(i).mean1-SC0(i).mean5,SC1(i).mean1-SC1(i).mean5); % group 1-5 vs permutation
%     SC1(i).p_diff_2v3 = SDL_p_permutation(SC0(i).mean2-SC0(i).mean3,SC1(i).mean2-SC1(i).mean3); % group 2-3 vs permutation
%     SC1(i).p_diff_2v4 = SDL_p_permutation(SC0(i).mean2-SC0(i).mean4,SC1(i).mean2-SC1(i).mean4); % group 2-4 vs permutation
%     SC1(i).p_diff_2v5 = SDL_p_permutation(SC0(i).mean2-SC0(i).mean5,SC1(i).mean2-SC1(i).mean5); % group 2-5 vs permutation
%     SC1(i).p_diff_3v4 = SDL_p_permutation(SC0(i).mean3-SC0(i).mean4,SC1(i).mean3-SC1(i).mean4); % group 3-4 vs permutation
%     SC1(i).p_diff_3v5 = SDL_p_permutation(SC0(i).mean3-SC0(i).mean5,SC1(i).mean3-SC1(i).mean5); % group 3-5 vs permutation
%     SC1(i).p_diff_4v5 = SDL_p_permutation(SC0(i).mean4-SC0(i).mean5,SC1(i).mean4-SC1(i).mean5); % group 4-5 vs permutation

   
    fprintf('Calculated: SC, CI & p-val for Top-%3d Areas\n',i);
end

%save
fdir = fullfile(SDL.out,SDL.data_type{1}); mkdir(fdir);
fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'SC0','SC1','-v7.3');
fprintf('Saved: Top-N analyses results saved in \n-->%s\n',fn);

% saved into a .mat fil, to be read by R
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn);
data = zeros(5001,147,8); % 5001(1 actual + 5000 random) x 147(network size: 2~148) x 20(very low/low/moderate/high/very high groups x raw/pos/neg/abs)
for j = 1:147 % per network size (2~148)
    data(1,j,1) = SC0(j+1).mean1raw; data(2:5001,j,1) = SC1(j+1).mean1raw; % mean1raw
    data(1,j,2) = SC0(j+1).mean2raw; data(2:5001,j,2) = SC1(j+1).mean2raw; % mean2raw
    data(1,j,3) = SC0(j+1).mean3raw; data(2:5001,j,3) = SC1(j+1).mean3raw; % mean3raw
    data(1,j,4) = SC0(j+1).mean4raw; data(2:5001,j,4) = SC1(j+1).mean4raw; % mean4raw
    data(1,j,5) = SC0(j+1).mean5raw; data(2:5001,j,5) = SC1(j+1).mean5raw; % mean5raw
    
    data(1,j,6) = SC0(j+1).mean1pos; data(2:5001,j,6) = SC1(j+1).mean1pos; % mean1pos
    data(1,j,7) = SC0(j+1).mean2pos; data(2:5001,j,7) = SC1(j+1).mean2pos; % mean2pos
    data(1,j,8) = SC0(j+1).mean3pos; data(2:5001,j,8) = SC1(j+1).mean3pos; % mean3pos
    data(1,j,9) = SC0(j+1).mean4pos; data(2:5001,j,9) = SC1(j+1).mean4pos; % mean4pos
    data(1,j,10)= SC0(j+1).mean5pos; data(2:5001,j,10)= SC1(j+1).mean5pos; % mean5pos
    
    data(1,j,11) = SC0(j+1).mean1neg; data(2:5001,j,11) = SC1(j+1).mean1neg; % mean1neg
    data(1,j,12) = SC0(j+1).mean2neg; data(2:5001,j,12) = SC1(j+1).mean2neg; % mean2neg
    data(1,j,13) = SC0(j+1).mean3neg; data(2:5001,j,13) = SC1(j+1).mean3neg; % mean3neg
    data(1,j,14) = SC0(j+1).mean4neg; data(2:5001,j,14) = SC1(j+1).mean4neg; % mean4neg
    data(1,j,15) = SC0(j+1).mean5neg; data(2:5001,j,15) = SC1(j+1).mean5neg; % mean5neg
    
    data(1,j,16) = SC0(j+1).mean1abs; data(2:5001,j,16) = SC1(j+1).mean1abs; % mean1abs
    data(1,j,17) = SC0(j+1).mean2abs; data(2:5001,j,17) = SC1(j+1).mean2abs; % mean2abs
    data(1,j,18) = SC0(j+1).mean3abs; data(2:5001,j,18) = SC1(j+1).mean3abs; % mean3abs
    data(1,j,19) = SC0(j+1).mean4abs; data(2:5001,j,19) = SC1(j+1).mean4abs; % mean4abs
    data(1,j,20) = SC0(j+1).mean5abs; data(2:5001,j,20) = SC1(j+1).mean5abs; % mean5abs
end
fn = fullfile(fdir,'Results_for_R.mat');
save(fn,'data');

% %% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT 
% fdir = fullfile(SDL.out,SDL.data_type{1}); 
% fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% % load(fn,'SC0','SC1');
% 
% X = 2:size(SC0,1);
% for i = 2:size(SC0,2)
%     
%     X1(i-1) = SC0(i).mean1; 
%     X2(i-1) = SC0(i).mean2; 
%     X3(i-1) = SC0(i).mean3; 
%     X4(i-1) = SC0(i).mean4; 
%     X5(i-1) = SC0(i).mean5; 
%     
%     p1(i-1) = SC1(i).p1; 
%     p2(i-1) = SC1(i).p2; 
%     p3(i-1) = SC1(i).p3; 
%     p4(i-1) = SC1(i).p4;
%     p5(i-1) = SC1(i).p5; 
%     
%     X1R(i-1) = mean(SC1(i).mean1); 
%     X2R(i-1) = mean(SC1(i).mean2); 
%     X3R(i-1) = mean(SC1(i).mean3);
%     X4R(i-1) = mean(SC1(i).mean4);
%     X5R(i-1) = mean(SC1(i).mean5);
%     
%     
%     
%     X1RA(i-1,:) = SC1(i).mean1; 
%     X2RA(i-1,:) = SC1(i).mean2; 
%     X3RA(i-1,:) = SC1(i).mean3; 
%     X4RA(i-1,:) = SC1(i).mean4; 
%     X5RA(i-1,:) = SC1(i).mean5; 
%     
% %     CI11L(i-1) = SC1(i).CI11(1); CI11R(i-1) = SC1(i).CI11(2); CI12L(i-1) = SC1(i).CI12(1); CI12R(i-1) = SC1(i).CI12(2);
% %     CI21L(i-1) = SC1(i).CI21(1); CI21R(i-1) = SC1(i).CI21(2); CI22L(i-1) = SC1(i).CI22(1); CI22R(i-1) = SC1(i).CI22(2);
% % 
% %     diff_11v21L(i-1) = SC1(i).diffCI_11v21(1); diff_11v21R(i-1) = SC1(i).diffCI_11v21(2);
% %     diff_12v22L(i-1) = SC1(i).diffCI_12v22(1); diff_12v22R(i-1) = SC1(i).diffCI_12v22(2);
% %     diff_11v12L(i-1) = SC1(i).diffCI_11v12(1); diff_11v12R(i-1) = SC1(i).diffCI_11v12(2);
% %     diff_21v22L(i-1) = SC1(i).diffCI_21v22(1); diff_21v22R(i-1) = SC1(i).diffCI_21v22(2);
% %     diff_main1L(i-1)  = SC1(i).diffCI_main1(1); diff_main1R(i-1) = SC1(i).diffCI_main1(2);
% %     diff_main2L(i-1)  = SC1(i).diffCI_main2(1); diff_main2R(i-1) = SC1(i).diffCI_main2(2);
% %     diff_interL(i-1)  = SC1(i).diffCI_inter(1); diff_interR(i-1) = SC1(i).diffCI_inter(2);
% %     
% %     
%     p_diff_1v2(i-1) = SC1(i).p_diff_1v2; 
%     p_diff_1v3(i-1) = SC1(i).p_diff_1v3; 
%     p_diff_1v4(i-1) = SC1(i).p_diff_1v4; 
%     p_diff_1v5(i-1) = SC1(i).p_diff_1v5; 
%     p_diff_2v3(i-1) = SC1(i).p_diff_2v3; 
%     p_diff_2v4(i-1) = SC1(i).p_diff_2v4; 
%     p_diff_2v5(i-1) = SC1(i).p_diff_2v5;
%     p_diff_3v4(i-1) = SC1(i).p_diff_3v4; 
%     p_diff_3v5(i-1) = SC1(i).p_diff_3v5;
%     p_diff_4v5(i-1) = SC1(i).p_diff_4v5;
% 
% end
% 
% 
% %% FDR corrected p values per top-area
% fprintf('\n\nFDR correction: group 1 vs random\n')
% [~, ~, ~, adj_p1] = fdr_bh(p1,0.05,'pdep','yes'); 
% find(adj_p1<=0.05)
% fprintf('\n\nFDR correction: group 2 vs random\n')
% [~, ~, ~, adj_p2] = fdr_bh(p2,0.05,'pdep','yes');
% find(adj_p2<=0.05)
% fprintf('\n\nFDR correction: group 3 vs random\n')
% [~, ~, ~, adj_p3] = fdr_bh(p3,0.05,'pdep','yes'); 
% find(adj_p3<=0.05)
% fprintf('\n\nFDR correction: group 4 vs random\n')
% [~, ~, ~, adj_p4] = fdr_bh(p4,0.05,'pdep','yes');
% find(adj_p4<=0.05)
% fprintf('\n\nFDR correction: group 5 vs random\n')
% [~, ~, ~, adj_p5] = fdr_bh(p5,0.05,'pdep','yes'); 
% find(adj_p5<=0.05)
% % 
% fprintf('\n\nFDR correction: group 1v2 vs random\n')
% [~, ~, ~, adj_p_diff_1v2] = fdr_bh(p_diff_1v2,0.05,'pdep','yes'); 
% find(adj_p_diff_1v2<=0.05)
% fprintf('\n\nFDR correction: group 1v3 vs random\n')
% [~, ~, ~, adj_p_diff_1v3] = fdr_bh(p_diff_1v3,0.05,'pdep','yes'); 
% find(adj_p_diff_1v3<=0.05)
% fprintf('\n\nFDR correction: group 1v4 vs random\n')
% [~, ~, ~, adj_p_diff_1v4] = fdr_bh(p_diff_1v4,0.05,'pdep','yes'); 
% find(adj_p_diff_1v4<=0.05)
% fprintf('\n\nFDR correction: group 1v5 vs random\n')
% [~, ~, ~, adj_p_diff_1v5] = fdr_bh(p_diff_1v5,0.05,'pdep','yes'); 
% find(adj_p_diff_1v5<=0.05)
% 
% fprintf('\n\nFDR correction: group 2v3 vs random\n')
% [~, ~, ~, adj_p_diff_2v3] = fdr_bh(p_diff_2v3,0.05,'pdep','yes'); 
% find(adj_p_diff_2v3<=0.05)
% fprintf('\n\nFDR correction: group 2v4 vs random\n')
% [~, ~, ~, adj_p_diff_2v4] = fdr_bh(p_diff_2v4,0.05,'pdep','yes'); 
% find(adj_p_diff_2v4<=0.05)
% fprintf('\n\nFDR correction: group 2v5 vs random\n')
% [~, ~, ~, adj_p_diff_2v5] = fdr_bh(p_diff_2v5,0.05,'pdep','yes'); 
% find(adj_p_diff_2v5<=0.05)
% 
% fprintf('\n\nFDR correction: group 3v4 vs random\n')
% [~, ~, ~, adj_p_diff_3v4] = fdr_bh(p_diff_3v4,0.05,'pdep','yes'); 
% find(adj_p_diff_3v4<=0.05)
% fprintf('\n\nFDR correction: group 3v5 vs random\n')
% [~, ~, ~, adj_p_diff_3v5] = fdr_bh(p_diff_3v5,0.05,'pdep','yes'); 
% find(adj_p_diff_3v5<=0.05)
% 
% fprintf('\n\nFDR correction: group 4v5 vs random\n')
% [~, ~, ~, adj_p_diff_4v5] = fdr_bh(p_diff_4v5,0.05,'pdep','yes'); 
% find(adj_p_diff_4v5<=0.05)
% 
% %% p value for area under curve
% p1a = SDL_p_permutation(mean(X1),mean(X1RA,1))
% p2a = SDL_p_permutation(mean(X2),mean(X2RA,1))
% p3a = SDL_p_permutation(mean(X3),mean(X3RA,1))
% p4a = SDL_p_permutation(mean(X4),mean(X4RA,1))
% p5a = SDL_p_permutation(mean(X5),mean(X5RA,1))
% 
% % 
% p_diff_1v2a = SDL_p_permutation(mean(X1)-mean(X2),mean(X1RA,1)-mean(X2RA,1))
% p_diff_1v3a = SDL_p_permutation(mean(X1)-mean(X3),mean(X1RA,1)-mean(X3RA,1))
% p_diff_1v4a = SDL_p_permutation(mean(X1)-mean(X4),mean(X1RA,1)-mean(X4RA,1))
% p_diff_1v5a = SDL_p_permutation(mean(X1)-mean(X5),mean(X1RA,1)-mean(X5RA,1))
% 
% p_diff_2v3a = SDL_p_permutation(mean(X2)-mean(X3),mean(X2RA,1)-mean(X3RA,1))
% p_diff_2v4a = SDL_p_permutation(mean(X2)-mean(X4),mean(X2RA,1)-mean(X4RA,1))
% p_diff_2v5a = SDL_p_permutation(mean(X2)-mean(X5),mean(X2RA,1)-mean(X5RA,1))
% 
% p_diff_3v4a = SDL_p_permutation(mean(X3)-mean(X4),mean(X3RA,1)-mean(X4RA,1))
% p_diff_3v5a = SDL_p_permutation(mean(X3)-mean(X5),mean(X3RA,1)-mean(X5RA,1))
% 
% p_diff_4v5a = SDL_p_permutation(mean(X4)-mean(X5),mean(X4RA,1)-mean(X5RA,1))
% 
% % multiple comparisons after FDR correction
% a = [p_diff_1v2a,p_diff_1v3a,p_diff_1v4a,p_diff_1v5a,p_diff_2v3a,p_diff_2v4a,p_diff_2v5a,p_diff_3v4a,p_diff_3v5a,p_diff_4v5a]
% [~, ~, ~, adj_p_diff_all] = fdr_bh(a,0.05,'pdep','yes')
% b = adj_p_diff_all;
% b = [1 b(1) b(2) b(3) b(4);...
%     b(1) 1 b(5) b(6) b(7);...
%     b(2),b(5),1,b(8),b(9);...
%     b(3),b(6),b(8),1,b(10);...
%     b(4),b(7),b(9),b(10),1]
% 
% 
% % get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
% % figure; % original vs permutation
% % L = 1:length(X1);
% % plot(L,(X1-X1R),...
% %     L,(X2-X2R),...
% %     L,(X3-X3R),...
% %     L,(X4-X4R),...
% %     L,(X5-X5R),'LineWidth', 2);
% % set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
% % a=prctile(T.CAPS,[20 40 60 80]);
% % a = [4 20 46 68]; % based on CT data
% % legend(['  0<=CAPS<',num2str(a(1))],...
% %     ['  ',num2str(a(1)),'<=CAPS<',num2str(a(2))],...
% %     [num2str(a(2)),'<=CAPS<',num2str(a(3))],...
% %     [num2str(a(3)),'<=CAPS<',num2str(a(4))],...
% %     [num2str(a(4)),'<=CAPS<120']);
% % xlabel('Top-N Regions','FontName', 'Arial','FontSize',16); set(gca,'XTick',0:10:150);
% % ylabel([SDL.data_type{1}(1:2),'-Based Structural Covariance'],'FontName', 'Arial','FontSize',16);
% % title('Real - Random','FontName', 'Arial','FontSize',16)
% % box off
% % 
% % 
% % 
% % b = sum([(X1-X1R)',...
% %     (X2-X2R)',...
% %     (X3-X3R)',...
% %     (X4-X4R)',...
% %     (X5-X5R)'])';
% % bar(0,b);
% % set(gca, 'FontName', 'Arial'); set(gca,'FontSize',16)
% % xlabel('PTSD Severity','FontName', 'Arial','FontSize',20); 
% % ylabel('AUC','FontName', 'Arial','FontSize',20);
% % xticklabels({' '});
% % box off
% 
% 
% % [x,y]=meshgrid(6,147);
% % z=[ (X11-X11R)-(X21-X21R);...
% %     (X12-X12R)-(X22-X22R);...
% %     (X13-X13R)-(X23-X23R);...
% %     (X14-X14R)-(X24-X24R);...
% %     (X15-X15R)-(X26-X26R);...
% %     (X16-X16R)-(X26-X26R)];
% % mesh(1:x,1:y,z)
% % colormap(jet)
% % title('Default colors for mesh BEFORE 2014b')
% % 
% % 
% % imagesc(z); colormap('hot'); colorbar;
% % xlabel('Top-N Regions'); ylabel('Age');title('Structural Covariance of (PTSD-CONT)Actual - (PTSD-CONT)Random')
% 
% 
% % savefig(fullfile(SDL.out,SDL.data_type{1},['Results_Curves_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));


%% End
end