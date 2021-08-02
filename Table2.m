function Table2(SDL)

% calculate values for Table 2, which is to show the PTSD vs random,
% non-PTSD vs random, and PTSD vs non-PTSD.
% (1) extract both global p-values and individual p-values from both CT-
% and SA-based outputs, and across three types of networks, i.e. atrophy,
% inflation and stable networks
% (2) FDR correction across 2 types of data (CT vs. SA) and 3 types of
% networks (atrophy, inflation vs. stable)

%% Parameters
SDL.data_type = {'CT';'SA'}; % 2 types of data
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
        fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{idata},'_corr.mat']);
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
        
        % mean of actual - random, 1x1
        m_global(inet,idata,1) = mean(X1-X1R);
        m_global(inet,idata,2) = mean(X2-X2R);

        
        % 95% CI (in fact 2.5~97.5% percentile), 1x2
        a = prctile(mean(X1'-X1RA),[2.5 97.5],2); err_global_L(inet,idata,1) = a(1); err_global_H(inet,idata,1) = a(2);
        a = prctile(mean(X2'-X2RA),[2.5 97.5],2); err_global_L(inet,idata,2) = a(1); err_global_H(inet,idata,2) = a(2);

        
        p_individual_PTSD(inet,idata,:) = p1;
        sign_individual_PTSD(inet,idata,:) = sign(X1-X1R); % sign of PTSD vs. random
        
        p_individual_nonPTSD(inet,idata,:) = p2;
        sign_individual_nonPTSD(inet,idata,:) = sign(X2-X2R); % sign of nonPTSD vs. random
        
        p_individual_PTSDvnonPTSD(inet,idata,:) = p_diff;
        sign_individual_PTSDvnonPTSD(inet,idata,:) = sign((X1-X1R)-(X2-X2R)); % sign of PTSD vs. nonPTSD
        
        
        p_global_PTSD(inet,idata)    = SDL_p_permutation(mean(X1),mean(X1RA,1));
        sign_global_PTSD(inet,idata) = sign(mean(X1)-mean(X1R)); % sign of PTSD vs. random
        
        p_global_nonPTSD(inet,idata)    = SDL_p_permutation(mean(X2),mean(X2RA,1));
        sign_global_nonPTSD(inet,idata) = sign(mean(X2)-mean(X2R)); % sign of nonPTSD vs. random
        
        if mean(X1)-mean(X2) > 0
            p_global_PTSDvnonPTSD(inet,idata) = SDL_p_permutation(mean(X1)-mean(X2),mean(X1RA,1)-mean(X2RA,1));
        else
            p_global_PTSDvnonPTSD(inet,idata) = SDL_p_permutation(mean(X2)-mean(X1),mean(X2RA,1)-mean(X1RA,1));
        end
        sign_global_PTSDvnonPTSD(inet,idata) = sign(mean(X1-X1R)-mean(X2-X2R)); % sign of PTSD vs. nonPTSD
        
    end
end


%% plot
figure
for i1 = 1:3
    for i2 = 1:2
        subplot(3,2,(i1-1)*2+i2);
        x = 1:2; 
        y = reshape(m_global(i1,i2,:),1,2); 
        errlow  = reshape(err_global_L(i1,i2,:),1,2); 
        errhigh = reshape(err_global_H(i1,i2,:),1,2);
        bar(x,y); 
        ylabel('Adj. Mean SC');
        set(gca,'xticklabel',{'PTSD','non-PTSD'}); 
        set(gca,'box','off'); 
        hold on; 
        er = errorbar(x,y,y-errlow,errhigh-y); 
        er.Color = [0 0 0]; 
        er.LineStyle = 'none'; hold off;
    end
end


%% Global p-values (FDR corrected across 2 data type and 3 network types)
[~, ~, ~, p_global_PTSD_adj] = fdr_bh(p_global_PTSD,0.05,'pdep','yes');
[~, ~, ~, p_global_nonPTSD_adj] = fdr_bh(p_global_nonPTSD,0.05,'pdep','yes');
[~, ~, ~, p_global_PTSDvnonPTSD_adj] = fdr_bh(p_global_PTSDvnonPTSD,0.05,'pdep','yes');

p_global_PTSD_adj = p_global_PTSD_adj .* sign_global_PTSD
p_global_nonPTSD_adj = p_global_nonPTSD_adj .* sign_global_nonPTSD
p_global_PTSDvnonPTSD_adj = p_global_PTSDvnonPTSD_adj .* sign_global_PTSDvnonPTSD


%% Individual p-values (FDR corrected across 2 data type, 3 network types, and 147 tests)
% corrected for 3 (types of networks) x 2 (types of data) x 147 (individual size)
[~, ~, ~, p_individual_PTSD_adj] = fdr_bh(p_individual_PTSD,0.05,'pdep','yes');
[~, ~, ~, p_individual_nonPTSD_adj] = fdr_bh(p_individual_nonPTSD,0.05,'pdep','yes');
[~, ~, ~, p_individual_PTSDvnonPTSD_adj] = fdr_bh(p_individual_PTSDvnonPTSD,0.05,'pdep','yes');

% a, the matrix to show both significant tests (P<0.05) and the sign of the tests (e.g. PTSD>random, or PTSD>nonPTSD)
% J, index of the significant tests
% J+1 because tests 1:147 corresponding to top-2 to -148 regions
% b is 1x1x107 for example, and b(:) is 107x1x1, consistent with J
fprintf('\n===Areas,individual_PTSD,FDR corrected\n'); a = (p_individual_PTSD_adj<=0.05).* sign_individual_PTSD;
for idata = 1:size(SDL.data_type,1),for inet = 1:size(SDL.net_type,1),fprintf('data=%s, net=%s, areas: ',SDL.data_type{idata},SDL.net_type{inet});
        J = find(a(inet,idata,:));b = a(inet,idata,J); fprintf('%d ',(J+1).*b(:));fprintf('\n');end,end

fprintf('\n===Areas,individual_nonPTSD,FDR corrected\n');a = (p_individual_nonPTSD_adj<=0.05).* sign_individual_nonPTSD;
for idata = 1:size(SDL.data_type,1),for inet = 1:size(SDL.net_type,1),fprintf('data=%s, net=%s, areas: ',SDL.data_type{idata},SDL.net_type{inet});
        J = find(a(inet,idata,:));b = a(inet,idata,J); fprintf('%d ',(J+1).*b(:));fprintf('\n');end,end

fprintf('\n===Areas,individual_PTSDvnonPTSD,FDR corrected\n');a = (p_individual_PTSDvnonPTSD_adj<=0.05).* sign_individual_PTSDvnonPTSD;
for idata = 1:size(SDL.data_type,1),for inet = 1:size(SDL.net_type,1),fprintf('data=%s, net=%s, areas: ',SDL.data_type{idata},SDL.net_type{inet});
        J = find(a(inet,idata,:));b = a(inet,idata,J); fprintf('%d ',(J+1).*b(:));fprintf('\n');end,end


% END
end