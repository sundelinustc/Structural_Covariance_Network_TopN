function SDL_power(SDL)

% 2020-02-18
% calculate sample size based on PTSD PGC data for Michael De Bellis's
% grant application (see below). We have CT and SA data for each age group
% 12, 13, 14, 15 and 16. We also know the PTSD vs Control group label. But
% we do not have the trauma info and drinking info
% ================
% Sorry to nag
% Any idea about power calculations for the grant, grant budgets folks are asking for a budget.
% How many would we need of youth with trauma (highly significant trauma like child abuse and neglect) in each age group age 12, 13, 14, 15, 16, I originally said 200 (40 in each group 50% gender mix).
%  
% The likihood of tradition to drinking will be 75% to 25% staying non drinkers by year 5, 60% to 35% in follow-up year 3.
%  
% Please advise, thanks Mike
%  
% Michael D. De Bellis, M.D., M.P.H.
% Professor of Psychiatry and Behavioral Sciences
% Director Healthy Childhood Brain Development and Developmental Traumatology Research Program, Department of Psychiatry and Behavioral Sciences
% Duke University Medical Center, Box 104360, Durham NC, 27710
% Phone: 919 812 3047
% Email: michael.debellis@duke.edu
% ================
% Planned analyses steps:
% (1) Load and clean data similiar to HBM manuscript Fig. 6C.
% (2) Use sampsizepwr.m to calculate sample size. The variables are the AUC (area under curve) of CT- or SA-based structural covariance in top-n regions. The variance comes from
%     the permuted data.
% (3) Plot power-sample size curves per age group

%% Load and clean data
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149}; % residules

idx11 = strcmp(T.Group,'PTSD') & (T.Age>=12 & T.Age<13); % group 11, PTSD & Age ~ 12
idx12 = strcmp(T.Group,'PTSD') & (T.Age>=13 & T.Age<14); % group 12, PTSD & Age ~ 13
idx13 = strcmp(T.Group,'PTSD') & (T.Age>=14 & T.Age<15); % group 13, PTSD & Age ~ 14
idx14 = strcmp(T.Group,'PTSD') & (T.Age>=15 & T.Age<16); % group 14, PTSD & Age ~ 15
idx15 = strcmp(T.Group,'PTSD') & (T.Age>=16 & T.Age<17); % group 15, PTSD & Age ~ 16

idx21 = strcmp(T.Group,'CONT') & (T.Age>=12 & T.Age<13); % group 21, CONT & Age ~ 12
idx22 = strcmp(T.Group,'CONT') & (T.Age>=13 & T.Age<14); % group 22, CONT & Age ~ 13
idx23 = strcmp(T.Group,'CONT') & (T.Age>=14 & T.Age<15); % group 23, CONT & Age ~ 14
idx24 = strcmp(T.Group,'CONT') & (T.Age>=15 & T.Age<16); % group 24, CONT & Age ~ 15
idx25 = strcmp(T.Group,'CONT') & (T.Age>=16 & T.Age<17); % group 25, CONT & Age ~ 16


[   sum(idx11),sum(idx21);...
    sum(idx12),sum(idx22);...
    sum(idx13),sum(idx23);...
    sum(idx14),sum(idx24);...
    sum(idx15),sum(idx25)]

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

% Structural Covariance Matrix -- Original     % SC of original PTSD & CONT data
S0 = SDL_SCA5(R,idx11,idx12,idx13,idx14,idx15,idx21,idx22,idx23,idx24,idx25,@corr,0); 

% loading
fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); % load SC0 and SC1 from existed results, i.e. CleanData, to avoid wasting time
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');
fprintf('Loaded: Top-N analyses midway results loaded from \n<--%s\n',fn);

% mean SC and std values
for i = 2:10:size(SC0,2) % per top-N
    % original data, values above diagonal
    node0 = SC0(i).node;
    SC0(i).Z11 = triu(S0.Y11(node0,node0),1); 
    SC0(i).Z12 = triu(S0.Y12(node0,node0),1); 
    SC0(i).Z13 = triu(S0.Y13(node0,node0),1);
    SC0(i).Z14 = triu(S0.Y14(node0,node0),1);
    SC0(i).Z15 = triu(S0.Y15(node0,node0),1);
    
    SC0(i).Z21 = triu(S0.Y21(node0,node0),1); 
    SC0(i).Z22 = triu(S0.Y22(node0,node0),1);
    SC0(i).Z23 = triu(S0.Y23(node0,node0),1);
    SC0(i).Z24 = triu(S0.Y24(node0,node0),1);
    SC0(i).Z25 = triu(S0.Y25(node0,node0),1);
    
    % non-zero values, Num of connections X 1
    SC0(i).Z11c = nonzeros(SC0(i).Z11); 
    SC0(i).Z12c = nonzeros(SC0(i).Z12); 
    SC0(i).Z13c = nonzeros(SC0(i).Z13);
    SC0(i).Z14c = nonzeros(SC0(i).Z14);
    SC0(i).Z15c = nonzeros(SC0(i).Z15);
    
    SC0(i).Z21c = nonzeros(SC0(i).Z21); 
    SC0(i).Z22c = nonzeros(SC0(i).Z22);% 
    SC0(i).Z23c = nonzeros(SC0(i).Z23);
    SC0(i).Z24c = nonzeros(SC0(i).Z24);
    SC0(i).Z25c = nonzeros(SC0(i).Z25);
    
    % mean of all connections, 1X1, original data
    SC0(i).mean11 = mean(SC0(i).Z11c); % mean of all connections, 1X1, original data, group 11
    SC0(i).mean12 = mean(SC0(i).Z12c); % mean of all connections, 1X1, original data, group 12
    SC0(i).mean13 = mean(SC0(i).Z13c);
    SC0(i).mean14 = mean(SC0(i).Z14c);
    SC0(i).mean15 = mean(SC0(i).Z15c);
    
    SC0(i).mean21 = mean(SC0(i).Z21c); % mean of all connections, 1X1, original data, group 21
    SC0(i).mean22 = mean(SC0(i).Z22c); % mean of all connections, 1X1, original data, group 22
    SC0(i).mean23 = mean(SC0(i).Z23c);
    SC0(i).mean24 = mean(SC0(i).Z24c);
    SC0(i).mean25 = mean(SC0(i).Z25c);

    
    % permuted data
    for j = 1:size(SC1(i).node,1) % per permutation (totally 5000 times)
        node1 = SC1(i).node(j,:);
        SC1(i).Z11(:,:,j) = triu(S0.Y11(node1,node1),1); 
        SC1(i).Z12(:,:,j) = triu(S0.Y12(node1,node1),1); 
        SC1(i).Z13(:,:,j) = triu(S0.Y13(node1,node1),1);
        SC1(i).Z14(:,:,j) = triu(S0.Y14(node1,node1),1);
        SC1(i).Z15(:,:,j) = triu(S0.Y15(node1,node1),1);
        
        SC1(i).Z21(:,:,j) = triu(S0.Y21(node1,node1),1); 
        SC1(i).Z22(:,:,j) = triu(S0.Y22(node1,node1),1);
        SC1(i).Z23(:,:,j) = triu(S0.Y23(node1,node1),1);
        SC1(i).Z24(:,:,j) = triu(S0.Y24(node1,node1),1);
        SC1(i).Z25(:,:,j) = triu(S0.Y25(node1,node1),1);
        
        SC1(i).Z11c(:,j) = nonzeros(SC1(i).Z11(:,:,j)); 
        SC1(i).Z12c(:,j) = nonzeros(SC1(i).Z12(:,:,j)); 
        SC1(i).Z13c(:,j) = nonzeros(SC1(i).Z13(:,:,j));
        SC1(i).Z14c(:,j) = nonzeros(SC1(i).Z14(:,:,j));
        SC1(i).Z15c(:,j) = nonzeros(SC1(i).Z15(:,:,j));
        
        SC1(i).Z21c(:,j) = nonzeros(SC1(i).Z21(:,:,j)); 
        SC1(i).Z22c(:,j) = nonzeros(SC1(i).Z22(:,:,j));% non-zero values, Num of connections X Num of permutations
        SC1(i).Z23c(:,j) = nonzeros(SC1(i).Z23(:,:,j));
        SC1(i).Z24c(:,j) = nonzeros(SC1(i).Z24(:,:,j));
        SC1(i).Z25c(:,j) = nonzeros(SC1(i).Z25(:,:,j));
    end
    SC1(i).mean11 = mean(SC1(i).Z11c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 11
    SC1(i).mean12 = mean(SC1(i).Z12c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 12
    SC1(i).mean13 = mean(SC1(i).Z13c,1);
    SC1(i).mean14 = mean(SC1(i).Z14c,1);
    SC1(i).mean15 = mean(SC1(i).Z15c,1);
    
    SC1(i).mean21 = mean(SC1(i).Z21c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 21
    SC1(i).mean22 = mean(SC1(i).Z22c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 22
    SC1(i).mean23 = mean(SC1(i).Z23c,1);
    SC1(i).mean24 = mean(SC1(i).Z24c,1);
    SC1(i).mean25 = mean(SC1(i).Z25c,1);
  
    fprintf('Calculated: SC, CI & p-val for Top-%3d Areas\n',i);
end
% %save
% fdir = fullfile(SDL.out,SDL.data_type{1}); mkdir(fdir);
% fn = fullfile(fdir,['PowerAnalyses_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% save(fn,'SC0','SC1','-v7.3');
% fprintf('Saved: Top-N analyses results saved in \n-->%s\n',fn);

%% Power analyses & plots
% fdir = fullfile(SDL.out,SDL.data_type{1}); 
% fn = fullfile(fdir,['PowerAnalyses_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% load(fn,'SC0','SC1');


X = 2:size(SC0,1);
j = 0;
for i = 2:10:148 % per calculation
    
    j = j + 1;
    
    X11(j) = SC0(i).mean11; 
    X12(j) = SC0(i).mean12; 
    X13(j) = SC0(i).mean13; 
    X14(j) = SC0(i).mean14; 
    X15(j) = SC0(i).mean15; 
    
    X21(j) = SC0(i).mean21; 
    X22(j) = SC0(i).mean22;
    X23(j) = SC0(i).mean23;
    X24(j) = SC0(i).mean24;
    X25(j) = SC0(i).mean25;

    
    
    
    X11RA(j,:) = SC1(i).mean11; 
    X12RA(j,:) = SC1(i).mean12; 
    X13RA(j,:) = SC1(i).mean13; 
    X14RA(j,:) = SC1(i).mean14; 
    X15RA(j,:) = SC1(i).mean15; 
    
    X21RA(j,:) = SC1(i).mean21; 
    X22RA(j,:) = SC1(i).mean22;
    X23RA(j,:) = SC1(i).mean23;
    X24RA(j,:) = SC1(i).mean24;
    X25RA(j,:) = SC1(i).mean25;

end


nn = 1:60; % sample size
nn1 = nn*1.15/0.9; % considering the +15% subjects for non-paremetric test, and 0.1 fail rate of scanning
ppow = 0.8; % power
figure;

% group 11-21, age = 12
p0.mean = mean(mean(X11RA-X21RA));% mean of difference under the null hypothesis
p0.std  =  std(mean(X11RA-X21RA));% std of  difference under the null hypothesis
p1 =      mean(X11-X21);% mean of difference under the alternative hypothesis
pwrout1 = sampsizepwr('t',[p0.mean p0.std],p1,[],nn); % power analysis group 11-21
nout1   = sampsizepwr('t',[p0.mean p0.std],p1,ppow)*1.15/0.9

% group 12-22, age = 13
p0.mean = mean(mean(X12RA-X22RA));% mean of difference under the null hypothesis
p0.std  =  std(mean(X12RA-X22RA));% std of  difference under the null hypothesis
p1 =      mean(X12-X22);% mean of difference under the alternative hypothesis
pwrout2 = sampsizepwr('t',[p0.mean p0.std],p1,[],nn); % power analysis group 12-22
nout2   = sampsizepwr('t',[p0.mean p0.std],p1,ppow)*1.15/0.9

% group 13-23, age = 14
p0.mean = mean(mean(X13RA-X23RA));% mean of difference under the null hypothesis
p0.std  =  std(mean(X13RA-X23RA));% std of  difference under the null hypothesis
p1 =      mean(X13-X23);% mean of difference under the alternative hypothesis
pwrout3 = sampsizepwr('t',[p0.mean p0.std],p1,[],nn); % power analysis group 13-23
nout3   = sampsizepwr('t',[p0.mean p0.std],p1,ppow)*1.15/0.9

% group 14-24, age = 15
p0.mean = mean(mean(X14RA-X24RA));% mean of difference under the null hypothesis
p0.std  =  std(mean(X14RA-X24RA));% std of  difference under the null hypothesis
p1 =      mean(X14-X24);% mean of difference under the alternative hypothesis
pwrout4 = sampsizepwr('t',[p0.mean p0.std],p1,[],nn); % power analysis group 14-24
nout4   = sampsizepwr('t',[p0.mean p0.std],p1,ppow)*1.15/0.9

% group 14-24, age = 15
p0.mean = mean(mean(X15RA-X25RA));% mean of difference under the null hypothesis
p0.std  =  std(mean(X15RA-X25RA));% std of  difference under the null hypothesis
p1 =      mean(X15-X25);% mean of difference under the alternative hypothesis
pwrout5 = sampsizepwr('t',[p0.mean p0.std],p1,[],nn); % power analysis group 15-25
nout5   = sampsizepwr('t',[p0.mean p0.std],p1,ppow)*1.15/0.9

plot(nn1,pwrout1,nn1,pwrout2,nn1,pwrout3);
% plot(nn,pwrout1,nn,pwrout2,nn,pwrout3,nn,pwrout4,nn,pwrout5);
title('Power versus Sample Size');xlabel('Sample Size');ylabel('Power')
legend('12 yrs','13 yrs','14 yrs', '15 yrs','16 yrs');


%% End
end