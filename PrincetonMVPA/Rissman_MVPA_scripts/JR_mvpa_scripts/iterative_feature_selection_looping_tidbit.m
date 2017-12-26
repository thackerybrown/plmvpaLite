num_vox = [25 50 100 250 500 750 1000 1500 2000 3000 4000 6000 8000 12000 16000 18000];

for j = [3 7 9 11 13 15]%1: length(num_vox)
    for i = 1:1      
        
        [subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_final_xval',['top_' num2str(num_vox(j)) '_voxels'],class_args);
        correct_vector = [];
        for a = 1:num_runs
            correct_vector = horzcat(correct_vector,results.iterations(a).perfmet.corrects);
        end
        overall_accuracy(i,j) = mean(correct_vector)
    end
end


class_args.train_funct_name = 'train_logreg';
class_args.test_funct_name = 'test_logreg';
class_args.penalty = 200;

class_args.train_funct_name = 'train_smlr';
class_args.test_funct_name = 'test_smlr';

class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';

num_vox = [25 50 100 250 500 750 1000 1500 2000 3000 4000 6000 8000 12000 16000 18000];

for j = 1: length(num_vox)
    for i = 1:1      
        
        [subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',['top_' num2str(num_vox(j)) '_voxels'],class_args);
        correct_vector = [];
        acts_diff_vector = [];
        for a = 1:5
            correct_vector = horzcat(correct_vector,results.iterations(a).perfmet.corrects);
            acts_diff_vector = horzcat(acts_diff_vector, results.iterations(a).acts(1,:)-results.iterations(a).acts(2,:));
        end
        
        [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
        abs_correct_sorted = correct_vector(abs_ind);
         num_trials = length(abs_correct_sorted);
             
        overall_accuracy(i,j) = mean(correct_vector)
        overall_top25pct(i,j) = mean(abs_correct_sorted(1:ceil(num_trials*.25)))
    end
end