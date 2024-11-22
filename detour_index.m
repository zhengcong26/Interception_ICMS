
function [detour_index_1,detour_index_2] = detour_index(classAll, t, nBootstraps)

for b = 1:nBootstraps
    
    for i = 1:size(classAll,1)
        N1 = size(classAll{i,1},1);
        idx1 = randi(N1, N1, 1);
        class_boot{i,1} = classAll{i,1}(idx1, :, :);
        
        N2 = size(classAll{i,2},1);
        idx2 = randi(N2, N1, 1);
        class_boot{i,2} = classAll{i,2}(idx2, :, :);
    end
    
    NSST_pca = cellfun(@(x) (squeeze(mean(x,1))),class_boot,'UniformOutput',0); % t * units
    
    NS_mat = cell2mat(NSST_pca(:,1));
    [coeff,score_CO,~,~,explained(:,b)] = pca(NS_mat);
    
    M = mean(NS_mat,1);
    INT_mat = cell2mat(NSST_pca(:,2));
    score_INT = (INT_mat-M)*coeff;
    
    score = [score_CO;score_INT];
    n = size(score,1)/4/2;
    nn = ones(1,8)*n;
    C = mat2cell(score,nn,size(score,2));
    
    ym_1 = C(1:4,:);
    ym_2 = C(5:8,:);
    
    %%
    if t(1)>-500
        t1 = find(t==-200);
    else
        t1 = find(t==-500);
    end
    t2 = find(t==100);
    
    pcn = 0;
    for PC = 1:length(explained(:,b))
        if sum(explained(1:PC,b))>90
            pcn = PC;
            break
        end
    end
    
%     pcn = 3;
    
    for k = 1:4
        mindit_1 = norm(ym_1{k,1}(t2,1:pcn)-ym_1{k,1}(t1,1:pcn));
        mindit_2 = norm(ym_2{k,1}(t2,1:pcn)-ym_2{k,1}(t1,1:pcn));
        
        sumdist_1=0;sumdist_2=0;
        for i = t1:t2-1
            sumdist_1 = sumdist_1+norm(ym_1{k,1}(i+1,1:pcn)-ym_1{k,1}(i,1:pcn));
            sumdist_2 = sumdist_2+norm(ym_2{k,1}(i+1,1:pcn)-ym_2{k,1}(i,1:pcn));
        end
        
        path_ratio_1(1,k)=sumdist_1/mindit_1;
        path_ratio_2(1,k)=sumdist_2/mindit_2;
    end
    
    detour_1(b,:) = path_ratio_1;
    detour_2(b,:) = path_ratio_2;
end

detour_index_1 = mean(detour_1,1);
detour_index_2 = mean(detour_2,1);

% mean(explained,2)



