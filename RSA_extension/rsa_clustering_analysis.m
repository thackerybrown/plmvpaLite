%% clustering
cm3 = corr(rmat_condensed);%generate a corrmat without NaNs


%% hierarchical clustering analysis

%create vector of distances between instances in the cm
distm = pdist(cm3,'correlation');% tell matlab metric is pearson r - this is critical; otherwise such algorithms assume values are some form of DISsimilarity
Z1 = linkage(distm,'average');%compute dendrogram, using average distance within clusters for agglomeration

%color code select original classes to help you evaluate clustering?
colorcodeZ1 = 1; %1 = yes, let's color code, 0 = no; you may turn this off if working out the color coding is unnecessary
if colorcodeZ1 == 1
    xz(1:size(names,2)) = {[0 0 0]};%create dendrogram condition colors using RGB code (default = [0 0 0], black)
    xz=xz';
    xz(logical(AA_intact)) = {[1 0 0]};
    xz(logical(EA_intact)) = {[0 1 0]};
    xz(logical(Scene_idx)) = {[0 0 1]};
    xz(logical(Obj_idx)) = {[1 0 1]};
    %h = cell2mat(userOptions.conditionColours)
    userOptions.conditionColours = cell2mat(xz);
end

% compute dendrogram. This example code labels each branch termination
% according to the name of that event in the 'names' model file you
% provided
figure;
subplot(1,1,1),[H_ignore, T_ignore, labelReordering] = dendrogram(Z1,0,'labels',names,'Orientation','left');%display dendrogram. 0 = show all items in correlation structure (default would limit to 30 clusters)

if colorcodeZ1 == 1
    color_t = xz(labelReordering);
    hold on;
    x = xlim(gca);
    for condition = 1:size(cm3,1)%size(squareRDM(cm4), 1)
        plot(x(1), condition, 'o', 'MarkerFaceColor', color_t{condition}, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    end%for:condition
    hold off;
end

% save plot
plot_savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '_hierarchicalclustering.png'];
saveas(gcf,plot_savename);

% ~~~~ that level of visualization can be too hard to make sense of

% Let's break it down with a couple of examples to simplify interpretation

% First, let's explore which classes belong to clusters at level __ in the dendrogram
clustlvl = 3;%define level of dendrogram to inspect. E.g., 3 means the level from the top where there are 3 clusters
cl_content{clustlvl} = cluster(Z1,'maxclust',clustlvl); %what are the items in cluster level 'cutoff'?
names_cl1 = names(cl_content{clustlvl}==1); % which classes are in cluster 1 at this level of the dendrogram?
names_cl2 = names(cl_content{clustlvl}==2); % which classes are in cluster 2 at this level of the dendrogram?
names_cl3 = names(cl_content{clustlvl}==3); % which classes are in cluster 3 at this level of the dendrogram?



%%% ok - let's visualize a simpler form to go along with that

% Let's try clustering just on a smaller subset of the data.
% Let's use the example from rsa_simple_2 using cm_2c and custlbls: that
% filtered correlation matrix is focused JUST on EA intact and AA intact -
% the question here, then, is can a clustering algorithm "discriminate"
% these two face types in our data?

%create vector of distances between instances in the cm
distm_2 = pdist(cm_2c,'correlation');% tell matlab metric is pearson r
Z2 = linkage(distm_2,'average');%compute dendrogram, using average distance within clusters for agglomeration
% compute dendrogram.
figure;
subplot(1,1,1),[H_ignore T_ignore labelReordering] = dendrogram(Z2,0,'labels',custlbls,'Orientation','left');%display dendrogram. 0 = show all items in correlation structure (default would limit to 30 clusters)
%color code select original classes to help you evaluate clustering?
colorcodeZ2 = 1; %1 = yes, let's color code, 0 = no; you may turn this off if working out the color coding is unnecessary
if colorcodeZ2 == 1
    %xz2(1:size(custlbls'),1) = {[0 0 0]};%create dendrogram condition colors using RGB code (default = [0 0 0], black)
    xz2(1:size(custlbls,2)) = {[0 0 0]};
    xz2 = xz2';
    xz2(logical(EA_intact(logical(testinds)))) = {[0 1 0]};%NOTE - before using our EA_intact idx here, we have to shorten it to just have 1s and 0s that agree with the size of our shrunken-down matrix; in this case, since our cm_2c drops values associated with AA, etc., we need to drop 0s in our EA_intact index that are associated with those same 0s. We can do this easily by filtering our EA_intact idx by the 'testinds' idx we created above for this example, which has 0s for every pattern which is dropped out of cm_2c and 1s for every pattern which is stil in cm_2c
    xz2(logical(Scene_idx(logical(testinds)))) = {[0 0 1]};
    userOptions.conditionColours2 = cell2mat(xz2);
    color_t2 = xz2(labelReordering);
    hold on;
    x = xlim(gca);
    for condition = 1:size(cm_2c,1)
        plot(x(1), condition, 'o', 'MarkerFaceColor', color_t2{condition}, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    end%for:condition
    hold off;
end

%%%% lastly - how do you know if your clustering fit the data well? 

% one way is to
% ensure the linkage heights it came up with resemble the raw pairwise
% distances in the original matrix.
% see, e.g., https://en.wikipedia.org/wiki/Cophenetic_correlation

c = cophenet(Z2,distm_2); % is c, the cophenetic correlation very high? ideally yes. if not, your linkage algorithm may not handle the data well and you might try something else (e.g., 'average' vs 'complete')

