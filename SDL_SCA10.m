function S = SDL_SCA10(R,idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx21,idx22,idx23,idx24,idx25,idx26,idx27,idx28,func,N)
% To calculate the mean SC and R-to-Z transformation
% Input
% R    --- Residules matrix after lme to regress out nusiance regressors
% idx11 ~ idx16 --- row index of group11 to group 16, i.e. PTSD & e.g. Age
%     in [0,20),[20,30),[30,40),[40,50),[50,60),[60,100)
% idx21 ~ idx26 --- row index of group21 to group 26, i.e. CONT & e.g. Age
%     in [0,20),[20,30),[30,40),[40,50),[50,60),[60,100)
% func --- matlab functions for calculation, i.e. @corr or @partialcorr
% N    --- number of permutations, no permutation if 0
% Output
% S    --- a structure containing
% S.M1 --- corr/partialcorr matrix of group1
% S.M2 --- corr/partialcorr matrix of group2
% S.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
% S.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2

S = [];
if N==0 % no permutation (for original data only)
    %% correlation or partial correlation
    S.M11 = func(R(idx11,:));
    S.M12 = func(R(idx12,:));
    S.M13 = func(R(idx13,:));
    S.M14 = func(R(idx14,:));
    S.M15 = func(R(idx15,:));
    S.M16 = func(R(idx16,:));
    S.M17 = func(R(idx17,:));
    S.M18 = func(R(idx18,:));
    
    S.M21 = func(R(idx21,:));
    S.M22 = func(R(idx22,:));
    S.M23 = func(R(idx23,:));
    S.M24 = func(R(idx24,:));
    S.M25 = func(R(idx25,:));
    S.M26 = func(R(idx26,:));
    S.M27 = func(R(idx27,:));
    S.M28 = func(R(idx28,:));
    %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
    S.Y11 = atanh(S.M11);
    S.Y12 = atanh(S.M12);
    S.Y13 = atanh(S.M13);
    S.Y14 = atanh(S.M14);
    S.Y15 = atanh(S.M15);
    S.Y16 = atanh(S.M16);
    S.Y17 = atanh(S.M17);
    S.Y18 = atanh(S.M18);
    
    S.Y21 = atanh(S.M21);
    S.Y22 = atanh(S.M22);
    S.Y23 = atanh(S.M23);
    S.Y24 = atanh(S.M24);
    S.Y25 = atanh(S.M25);
    S.Y26 = atanh(S.M26);
    S.Y27 = atanh(S.M27);
    S.Y28 = atanh(S.M28);
    
%     get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
%     figure;
%     subplot(2,4,1);histogram(S.M11);  xlabel('Coeff.'); ylabel('Count');title('Group = PTSD &1');
%     subplot(2,4,2);histogram(S.M12);  xlabel('Coeff.'); ylabel('Count');title('Group = PTSD &2');
%     subplot(2,4,3);histogram(S.M21);  xlabel('Coeff.'); ylabel('Count');title('Group = CONT &1');
%     subplot(2,4,4);histogram(S.M22);  xlabel('Coeff.'); ylabel('Count');title('Group = CONT &2');
%     subplot(2,4,5);histogram(S.Y11);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,6);histogram(S.Y12);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,7);histogram(S.Y21);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,8);histogram(S.Y22);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
    fprintf('Completed: SC matrix and R-to-Z transformation, N=0\n\n');
else % permutation
%     idx0 = [idx11;idx12;idx21;idx22]; % shuffling group labels 
%     for i = 1:N
%         %% Shuffle labels
%         idx = idx0(randperm(length(idx0))); % shuffle the labels of two groups
%         %% correlation or partial correlation
%         S.M11(:,:,i) = func(R(idx(idx11),:));
%         S.M12(:,:,i) = func(R(idx(idx12),:));
%         S.M13(:,:,i) = func(R(idx(idx13),:));
%         S.M14(:,:,i) = func(R(idx(idx14),:));
%         S.M15(:,:,i) = func(R(idx(idx15),:));
%         S.M16(:,:,i) = func(R(idx(idx16),:));
%         
%         S.M21(:,:,i) = func(R(idx(idx21),:));
%         S.M22(:,:,i) = func(R(idx(idx22),:));
%         S.M23(:,:,i) = func(R(idx(idx23),:));
%         S.M24(:,:,i) = func(R(idx(idx24),:));
%         S.M25(:,:,i) = func(R(idx(idx25),:));
%         S.M26(:,:,i) = func(R(idx(idx26),:));
%         %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
%         S.Y11 = atanh(S.M11);
%         S.Y12 = atanh(S.M12);
%         S.Y13 = atanh(S.M13);
%         S.Y14 = atanh(S.M14);
%         S.Y15 = atanh(S.M15);
%         S.Y16 = atanh(S.M16);
%         
%         S.Y21 = atanh(S.M21);
%         S.Y22 = atanh(S.M22);
%         S.Y23 = atanh(S.M23);
%         S.Y24 = atanh(S.M24);
%         S.Y25 = atanh(S.M25);
%         S.Y26 = atanh(S.M26);
%         fprintf('Completed: SC matrix and R-to-Z transformation, N=%4d\n',i);
%     end
end
end