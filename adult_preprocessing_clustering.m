%% import table
%clear
adult = readtable('adult.csv');
adult.Properties.VariableNames = {'age', 'workclass', 'fnlwgt',...
    'education', 'education_num', 'marital_status', 'occupation', ...
    'relationship', 'race', 'sex', 'capital_gain', 'capital_loss', ...
    'hours_per_week', 'native_country', 'result'};
adult = removevars(adult,{'fnlwgt'}); %not important, it indicates how widespread is this kind of person in the census
adult = removevars(adult,{'education_num'}); %not important either, education is already encoded

adult = standardizeMissing(adult, {'?'}); %find and replace ? with missing value

cate_feats = {'workclass', 'education', 'marital_status', 'occupation', ...
    'relationship', 'race', 'sex','native_country','result'};
num_feats = {'age','capital_gain',...
    'capital_loss','hours_per_week'};
adult_rencoded = table();

for j=1:length(num_feats)
    adult_rencoded(:,num_feats{j}) = array2table(...
        table2array(adult(:,num_feats{j})),'VariableNames',{num_feats{j}});
end

for j=1:length(cate_feats)
    adult_rencoded(:,cate_feats{j}) = array2table(categorical(...
        table2array(adult(:,cate_feats{j}))),'VariableNames',{cate_feats{j}});
end

adult_rencoded.result = (grp2idx(adult_rencoded.result) - 1) > 0 ;

findnans=ismissing(adult_rencoded);
find_rows=sum(findnans,2)>0; %find rows with missing values
adult_full = adult_rencoded(~find_rows,:);
cate_feats = cate_feats(1:8);
%%
[grouped_clusters, lone_clusters] = cutdown_similars(adult_full,cate_feats, num_feats, 0.99999999999);
save grouped_clusters.mat grouped_clusters lone_clusters

%%
% TODO ripeti la procedura con un numero di massimi pari a 8 per le
% eccezioni già trovate
% fermati fino a che non consideri le eccezioni ripulite
[lone_grouped_clusters, loner_clusters] = cutdown_similars(lone_clusters,cate_feats, num_feats, 0.9999);
%%
function [grouped_clusters, lone_clusters] = cutdown_similars(adult_full, cate_feats, num_feats, var)
    % onehot encode the categorical features
    adult_new = adult_full;
    adult_cat = adult_new(:,cate_feats);
    adult_cat_ohe = table();
    for i=1:8
        temp = onehotencode(adult_cat(:,cate_feats{i}));
        adult_cat_ohe = [adult_cat_ohe,temp];
    end
    % make clusters out of the numerical features
    max_found = length(cate_feats) + 1;
    grouped_clusters = table();
    iter = 0;
    
    while max_found >= 8
    %while iter<10
        tic
        % numeric features
        adult_numeric = adult_new(:,num_feats);
        
        % score is the actual indexes to make fcm
        [~,score,~,~,explained] = pca(table2array(adult_numeric));
        
        prop = 0; variance = 0;
        
        while variance<var
            prop = prop+1;
            variance = sum(explained(1:prop));
        end
        
        disp(strcat('Number of values: ', num2str(prop)))
        [~, U] = fcm(score,2);
        [~,imaxU] = max(U);
        clust_names = cell(1,prop);
        
        for i=1:(prop+1)
            clust_names(i) = {strcat('c_',num2str(i))};
        end
        
        clusters = categorical(imaxU',1:(prop+1),clust_names);
        cluster_ohe = onehotencode(clusters,2);
        disp(size(cluster_ohe))
        disp(size(adult_cat_ohe))
        adult_full_ohe = [cluster_ohe, table2array(adult_cat_ohe)];
        
        disp('calculating equality matrix')
        
        eq_matrix = adult_full_ohe*adult_full_ohe';
        
        disp('removing diagonal elements')
        eq_mat = eq_matrix.*~eye(size(eq_matrix));
        
        max_found = max(eq_mat,[],'all'); % this changes often
        max_eq = eq_mat==max_found;
        identic = sum(max_eq);
        
        disp(strcat('Max similarity found: ',num2str(max_found)))
        [~,I] = max(identic);
        
        % change the values 
        newgroup_index = max_eq(:,I);
        newgroup_index(I) = true;
        newgroup = adult_new(newgroup_index,:);
        newgroup(:,'group_index') = array2table(iter*ones(height(newgroup),1),'VariableNames',{'group_index'});
        disp(strcat('Number of rows added to the cleaned dataset: ',num2str(height(newgroup))))
        
        grouped_clusters = [grouped_clusters;newgroup];
        
        disp(strcat('Number of rows of the cleaned dataset: ',num2str(height(grouped_clusters))))
        
        adult_new=adult_new(~newgroup_index,:); 
        
        disp(strcat('Number of rows of the cleaned dataset: ',num2str(height(adult_new))))
        disp(height(grouped_clusters)-height(newgroup))
        adult_cat_ohe = adult_cat_ohe(~newgroup_index,:);
        iter = iter + 1;
        toc
    end
    lone_clusters = adult_new;
end