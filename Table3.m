function Table3(SDL)

% calculate values for Table 3, which is to show the PTSD severity.
% (1) extract both global p-values and individual p-values from both CT-
% and SA-based outputs, and across three types of networks, i.e. atrophy,
% inflation and stable networks
% (2) FDR correction across 2 types of data (CT vs. SA) and 3 types of
% networks (atrophy, inflation vs. stable) and 10 between-group comparisons

%% Parameters
SDL.data_type = {'CT_PTSDsev10';'SA_PTSDsev10'}; % 2 types of data
SDL.net_type  = {'Reduction';'Increase';'Nodiff'}; % 3 types of networks

for idata = 1:size(SDL.data_type,1) % per data type
    for inet = 1:size(SDL.net_type,1) % per network type
        
        %% Loading information
        % Residuals
        SDL.out = fullfile(SDL.path,['Outputs_ComBat_',SDL.net_type{inet}]);
        fdir = fullfile(SDL.out,SDL.data_type{idata});
        fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{idata},'_corr_.mat']);
        load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
        R = T{:,2:149}; % residules
        idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups
        
        % Effect size
        fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{idata},'_corr_.mat']);
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
        S0 = SDL_SCA(R,idx1,idx2,@corr,0);     % SC of original PTSD & CONT data
        % S0.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
        % S0.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2
        
        %% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT
        fdir = fullfile(SDL.out,SDL.data_type{idata});
        fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{idata},'_corr.mat']);
        load(fn,'SC0','SC1');
        
        X = 2:size(SC0,1);
        for i = 2:size(SC0,2)
            
            X1(i-1) = SC0(i).mean1;
            X2(i-1) = SC0(i).mean2;
            X3(i-1) = SC0(i).mean3;
            X4(i-1) = SC0(i).mean4;
            X5(i-1) = SC0(i).mean5;
            
            p1(i-1) = SC1(i).p1;
            p2(i-1) = SC1(i).p2;
            p3(i-1) = SC1(i).p3;
            p4(i-1) = SC1(i).p4;
            p5(i-1) = SC1(i).p5;
            
            X1R(i-1) = mean(SC1(i).mean1);
            X2R(i-1) = mean(SC1(i).mean2);
            X3R(i-1) = mean(SC1(i).mean3);
            X4R(i-1) = mean(SC1(i).mean4);
            X5R(i-1) = mean(SC1(i).mean5);
            
            
            
            X1RA(i-1,:) = SC1(i).mean1;
            X2RA(i-1,:) = SC1(i).mean2;
            X3RA(i-1,:) = SC1(i).mean3;
            X4RA(i-1,:) = SC1(i).mean4;
            X5RA(i-1,:) = SC1(i).mean5;
            
            p_diff_1v2(i-1) = SC1(i).p_diff_1v2;
            p_diff_1v3(i-1) = SC1(i).p_diff_1v3;
            p_diff_1v4(i-1) = SC1(i).p_diff_1v4;
            p_diff_1v5(i-1) = SC1(i).p_diff_1v5;
            p_diff_2v3(i-1) = SC1(i).p_diff_2v3;
            p_diff_2v4(i-1) = SC1(i).p_diff_2v4;
            p_diff_2v5(i-1) = SC1(i).p_diff_2v5;
            p_diff_3v4(i-1) = SC1(i).p_diff_3v4;
            p_diff_3v5(i-1) = SC1(i).p_diff_3v5;
            p_diff_4v5(i-1) = SC1(i).p_diff_4v5;
            
        end
        
        %% individual values
        % mean of actual - random, 3(nets)x2(data)x5(groups)x147(comparisons)
        m_individual_diff(inet,idata,1,:) = X1-X1R; 
        m_individual_diff(inet,idata,2,:) = X2-X2R;
        m_individual_diff(inet,idata,3,:) = X3-X3R;
        m_individual_diff(inet,idata,4,:) = X4-X4R;
        m_individual_diff(inet,idata,5,:) = X5-X5R;
        
        % p-values, random, 3(nets)x2(data)x10(between-group contrasts)x147(comparisons)
        p_individual_diff(inet,idata,1,:) = p_diff_1v2; 
        p_individual_diff(inet,idata,2,:) = p_diff_1v3; 
        p_individual_diff(inet,idata,3,:) = p_diff_1v4; 
        p_individual_diff(inet,idata,4,:) = p_diff_1v5; 
        p_individual_diff(inet,idata,5,:) = p_diff_2v3; 
        p_individual_diff(inet,idata,6,:) = p_diff_2v4; 
        p_individual_diff(inet,idata,7,:) = p_diff_2v5; 
        p_individual_diff(inet,idata,8,:) = p_diff_3v4; 
        p_individual_diff(inet,idata,9,:) = p_diff_3v5; 
        p_individual_diff(inet,idata,10,:)= p_diff_4v5; 
        
        % sign, random, 3(nets)x2(data)x10(between-group contrasts)x147(comparisons)
        sign_individual_diff(inet,idata,1,:) = sign((X1-X1R)-(X2-X2R)); 
        sign_individual_diff(inet,idata,2,:) = sign((X1-X1R)-(X3-X3R));
        sign_individual_diff(inet,idata,3,:) = sign((X1-X1R)-(X4-X4R));
        sign_individual_diff(inet,idata,4,:) = sign((X1-X1R)-(X5-X5R));
        sign_individual_diff(inet,idata,5,:) = sign((X2-X2R)-(X3-X3R));
        sign_individual_diff(inet,idata,6,:) = sign((X2-X2R)-(X4-X4R));
        sign_individual_diff(inet,idata,7,:) = sign((X2-X2R)-(X5-X5R));
        sign_individual_diff(inet,idata,8,:) = sign((X3-X3R)-(X4-X4R));
        sign_individual_diff(inet,idata,9,:) = sign((X3-X3R)-(X5-X5R));
        sign_individual_diff(inet,idata,10,:)= sign((X4-X4R)-(X5-X5R));
        
        
        %%  global values
        % mean of actual - random, 1x1
        m_global(inet,idata,1) = mean(X1-X1R);
        m_global(inet,idata,2) = mean(X2-X2R);
        m_global(inet,idata,3) = mean(X3-X3R);
        m_global(inet,idata,4) = mean(X4-X4R);
        m_global(inet,idata,5) = mean(X5-X5R);
        
        % 95% CI (in fact 2.5~97.5% percentile), 1x2
        CI_global_1(inet,idata,:) = prctile(mean(X1'-X1RA),[2.5 97.5],2); err_global_L(inet,idata,1) =  CI_global_1(inet,idata,1); err_global_H(inet,idata,1) =  CI_global_1(inet,idata,2);
        CI_global_2(inet,idata,:) = prctile(mean(X2'-X2RA),[2.5 97.5],2); err_global_L(inet,idata,2) =  CI_global_2(inet,idata,1); err_global_H(inet,idata,2) =  CI_global_2(inet,idata,2);
        CI_global_3(inet,idata,:) = prctile(mean(X3'-X3RA),[2.5 97.5],2); err_global_L(inet,idata,3) =  CI_global_3(inet,idata,1); err_global_H(inet,idata,3) =  CI_global_3(inet,idata,2);
        CI_global_4(inet,idata,:) = prctile(mean(X4'-X4RA),[2.5 97.5],2); err_global_L(inet,idata,4) =  CI_global_4(inet,idata,1); err_global_H(inet,idata,4) =  CI_global_4(inet,idata,2);
        CI_global_5(inet,idata,:) = prctile(mean(X5'-X5RA),[2.5 97.5],2); err_global_L(inet,idata,5) =  CI_global_5(inet,idata,1); err_global_H(inet,idata,5) =  CI_global_5(inet,idata,2);
        
        
        % p values, 1x1
        if mean(X1-X2)>0, p_global_diff_1v2(inet,idata)=SDL_p_permutation(mean(X1-X2),mean(X1RA-X2RA,1)); else,p_global_diff_1v2(inet,idata)=SDL_p_permutation(mean(X2-X1),mean(X2RA-X1RA,1));end
        if mean(X1-X3)>0, p_global_diff_1v3(inet,idata)=SDL_p_permutation(mean(X1-X3),mean(X1RA-X3RA,1)); else,p_global_diff_1v3(inet,idata)=SDL_p_permutation(mean(X3-X1),mean(X3RA-X1RA,1));end
        if mean(X1-X4)>0, p_global_diff_1v4(inet,idata)=SDL_p_permutation(mean(X1-X4),mean(X1RA-X4RA,1)); else,p_global_diff_1v4(inet,idata)=SDL_p_permutation(mean(X4-X1),mean(X4RA-X1RA,1));end
        if mean(X1-X5)>0, p_global_diff_1v5(inet,idata)=SDL_p_permutation(mean(X1-X5),mean(X1RA-X5RA,1)); else,p_global_diff_1v5(inet,idata)=SDL_p_permutation(mean(X5-X1),mean(X5RA-X1RA,1));end
        if mean(X2-X3)>0, p_global_diff_2v3(inet,idata)=SDL_p_permutation(mean(X2-X3),mean(X2RA-X3RA,1)); else,p_global_diff_2v3(inet,idata)=SDL_p_permutation(mean(X3-X2),mean(X3RA-X2RA,1));end
        if mean(X2-X4)>0, p_global_diff_2v4(inet,idata)=SDL_p_permutation(mean(X2-X4),mean(X2RA-X4RA,1)); else,p_global_diff_2v4(inet,idata)=SDL_p_permutation(mean(X4-X2),mean(X4RA-X2RA,1));end
        if mean(X2-X5)>0, p_global_diff_2v5(inet,idata)=SDL_p_permutation(mean(X2-X5),mean(X2RA-X5RA,1)); else,p_global_diff_2v5(inet,idata)=SDL_p_permutation(mean(X5-X2),mean(X5RA-X2RA,1));end
        if mean(X3-X4)>0, p_global_diff_3v4(inet,idata)=SDL_p_permutation(mean(X3-X4),mean(X3RA-X4RA,1)); else,p_global_diff_3v4(inet,idata)=SDL_p_permutation(mean(X4-X3),mean(X4RA-X3RA,1));end
        if mean(X3-X5)>0, p_global_diff_3v5(inet,idata)=SDL_p_permutation(mean(X3-X5),mean(X3RA-X5RA,1)); else,p_global_diff_3v5(inet,idata)=SDL_p_permutation(mean(X5-X3),mean(X5RA-X3RA,1));end
        if mean(X4-X5)>0, p_global_diff_4v5(inet,idata)=SDL_p_permutation(mean(X4-X5),mean(X4RA-X5RA,1)); else,p_global_diff_4v5(inet,idata)=SDL_p_permutation(mean(X5-X4),mean(X5RA-X4RA,1));end
        
        % sign
        sign_global_diff_1v2(inet,idata) = sign(mean((X1-X1R)-(X2-X2R)));
        sign_global_diff_1v3(inet,idata) = sign(mean((X1-X1R)-(X3-X3R)));
        sign_global_diff_1v4(inet,idata) = sign(mean((X1-X1R)-(X4-X4R)));
        sign_global_diff_1v5(inet,idata) = sign(mean((X1-X1R)-(X5-X5R)));
        sign_global_diff_2v3(inet,idata) = sign(mean((X2-X2R)-(X3-X3R)));
        sign_global_diff_2v4(inet,idata) = sign(mean((X2-X2R)-(X4-X4R)));
        sign_global_diff_2v5(inet,idata) = sign(mean((X2-X2R)-(X5-X5R)));
        sign_global_diff_3v4(inet,idata) = sign(mean((X3-X3R)-(X4-X4R)));
        sign_global_diff_3v5(inet,idata) = sign(mean((X3-X3R)-(X5-X5R)));
        sign_global_diff_4v5(inet,idata) = sign(mean((X4-X4R)-(X5-X5R)));
        
        
    end
end


%% plot
figure
for i1 = 1:3
    for i2 = 1:2
        subplot(3,2,(i1-1)*2+i2);
        x = 1:5; y = reshape(m_global(i1,i2,:),1,5); errlow = reshape(err_global_L(i1,i2,:),1,5); errhigh = reshape(err_global_H(i1,i2,:),1,5);
        bar(x,y); set(gca,'xticklabel',{'g1','g2','g3','g4','g5'}); set(gca,'box','off'); hold on; er = errorbar(x,y,y-errlow,errhigh-y); er.Color = [0 0 0]; er.LineStyle = 'none'; hold off;
    end
end


%% Global p-values (FDR corrected across 2 data type, 3 network types and 10 comparisons)
p_global_diff_all(:,:,1) = p_global_diff_1v2; sign_global_diff_all(:,:,1) =  sign_global_diff_1v2;
p_global_diff_all(:,:,2) = p_global_diff_1v3; sign_global_diff_all(:,:,2) =  sign_global_diff_1v3;
p_global_diff_all(:,:,3) = p_global_diff_1v4; sign_global_diff_all(:,:,3) =  sign_global_diff_1v4;
p_global_diff_all(:,:,4) = p_global_diff_1v5; sign_global_diff_all(:,:,4) =  sign_global_diff_1v5;
p_global_diff_all(:,:,5) = p_global_diff_2v3; sign_global_diff_all(:,:,5) =  sign_global_diff_2v3;
p_global_diff_all(:,:,6) = p_global_diff_2v4; sign_global_diff_all(:,:,6) =  sign_global_diff_2v4;
p_global_diff_all(:,:,7) = p_global_diff_2v5; sign_global_diff_all(:,:,7) =  sign_global_diff_2v5;
p_global_diff_all(:,:,8) = p_global_diff_3v4; sign_global_diff_all(:,:,8) =  sign_global_diff_3v4;
p_global_diff_all(:,:,9) = p_global_diff_3v5; sign_global_diff_all(:,:,9) =  sign_global_diff_3v5;
p_global_diff_all(:,:,10)= p_global_diff_4v5; sign_global_diff_all(:,:,10)=  sign_global_diff_4v5;
[~, ~, ~, p_global_diff_all_adj] = fdr_bh(p_global_diff_all,0.05,'pdep','yes');
p_global_diff_all_adj = p_global_diff_all_adj.*sign_global_diff_all; % p-values with sign
for idata = 1:size(SDL.data_type,1) % per data type
    for inet = 1:size(SDL.net_type,1) % per network type
        p_global_diff(inet,idata,1:5,1:5) = 1; % initial values, NS
        p_global_diff(inet,idata,1,2) = p_global_diff_all_adj(inet,idata,1); p_global_diff(inet,idata,2,1) = p_global_diff_all_adj(inet,idata,1);
        p_global_diff(inet,idata,1,3) = p_global_diff_all_adj(inet,idata,2); p_global_diff(inet,idata,3,1) = p_global_diff_all_adj(inet,idata,2);
        p_global_diff(inet,idata,1,4) = p_global_diff_all_adj(inet,idata,3); p_global_diff(inet,idata,4,1) = p_global_diff_all_adj(inet,idata,3);
        p_global_diff(inet,idata,1,5) = p_global_diff_all_adj(inet,idata,4); p_global_diff(inet,idata,5,1) = p_global_diff_all_adj(inet,idata,4);
        p_global_diff(inet,idata,2,3) = p_global_diff_all_adj(inet,idata,5); p_global_diff(inet,idata,3,2) = p_global_diff_all_adj(inet,idata,5);
        p_global_diff(inet,idata,2,4) = p_global_diff_all_adj(inet,idata,6); p_global_diff(inet,idata,4,2) = p_global_diff_all_adj(inet,idata,6);
        p_global_diff(inet,idata,2,5) = p_global_diff_all_adj(inet,idata,7); p_global_diff(inet,idata,5,2) = p_global_diff_all_adj(inet,idata,7);
        p_global_diff(inet,idata,3,4) = p_global_diff_all_adj(inet,idata,8); p_global_diff(inet,idata,4,3) = p_global_diff_all_adj(inet,idata,8);
        p_global_diff(inet,idata,3,5) = p_global_diff_all_adj(inet,idata,9); p_global_diff(inet,idata,5,3) = p_global_diff_all_adj(inet,idata,9);
        p_global_diff(inet,idata,4,5) =p_global_diff_all_adj(inet,idata,10); p_global_diff(inet,idata,5,4) = p_global_diff_all_adj(inet,idata,10);
    end
end


%% Individual p-values (FDR corrected across 2 data type, 3 network types, and 147 tests)
% corrected for 3 (types of networks) x 2 (types of data) x 147 (individual size)
[~, ~, ~, p_individual_diff_adj] = fdr_bh(p_individual_diff,0.05,'pdep','yes');
% a, the matrix to show both significant tests (P<0.05) and the sign of the tests (e.g. PTSD>random, or PTSD>nonPTSD)
% J, index of the significant tests
% J+1 because tests 1:147 corresponding to top-2 to -148 regions
% b is 1x1x107 for example, and b(:) is 107x1x1, consistent with J

fprintf('\n===Areas,individual_diff,FDR corrected\n');a = (p_individual_diff_adj<=0.05).* sign_individual_diff;
for idata = 1:size(SDL.data_type,1),for inet = 1:size(SDL.net_type,1),for icon = 1:10
        fprintf('data=%s, net=%s, con=%d, areas: ',SDL.data_type{idata},SDL.net_type{inet},icon);
        J = find(a(inet,idata,icon,:));b = a(inet,idata,icon,J); fprintf('%d ',(J+1).*b(:));fprintf('\n');
    end,end,end


% END
end