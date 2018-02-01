function CustomIE_TIB(f_path, str2blowup, leftover)                            
% input the file's location name as a string, the cell of condition names
% in the order you you want output (to expand), and the cell of unwanted
% conditions in the order to keep track (will not have onsets or durations
% expanded.
%example: CustomIE_TIB('C:\Users\user\Documents\Work\circmaze\modelfiles\onsets_cm25_allruns.mat', {'cues1', 'cues2', 'cues3', 'cues4', 'cues5', 'goal1', 'goal2', 'goal3', 'goal4', 'goal5'},{'plan1', 'plan2', 'plan3', 'plan4', 'plan5', 'nav1', 'nav2', 'nav3', 'nav4', 'nav5', 'goalhold1', 'goalhold2', 'goalhold3', 'goalhold4', 'goalhold5', 'arrows'})

 load(f_path);  % loads file                                               % 
    dur = []; % makes empty matrices for later concatenation
    on = [];
    
    for j = 1:length(str2blowup) %this returns the total number of variables per cell in the cell array in a similar matrices
        dur = [dur length(durations{j})];
        on = [on length(onsets{j})];
    end
    
    dur2 = cell(1,0); % again, a cell array made to make concatenation easier
    on2 = cell(1,0);
    nam2 = cell(1,0);
    
for i = 1:length(str2blowup)                                                % this loop goes through the cell of wanted conditions, matches them, and expands duration,
    ind = cellfun(@(c) strcmp(c,str2blowup{i}),names);                            % onsets, and names in the order you fed the function      
    inds2blowup = find(ind==1); % indices where names match in previous format
    dur2 = cat(2,dur2,num2cell(durations{inds2blowup}));                    
    on2 = cat(2,on2,num2cell(onsets{inds2blowup}));                        % concatenation of wanted durations and onsets for particular condition
    nam = cell(1,length(onsets{inds2blowup}));                          % preallocates a cell for the name expansion
        for j = 1:length(onsets{inds2blowup})                           % this loop assigns the correct name and number of names
            nam{j} = names{ind};                                                             % to the condition durations/onsets indices
        end
    nam2 = cat(2,nam2,nam);                                                % tacks on the names just made to overall names list
%     dur = cell(1,length(onsets{inds2blowup}));                          % preallocates a cell for the name expansion
%         for j = 1:length(onsets{inds2blowup})                           % this loop assigns the correct name and number of names
%             dur{j} = durations{ind};                                                             % to the condition durations/onsets indices
%         end
%     dur2 = cat(2,dur2,dur);        
end
dur3 = cell(1,length(leftover)); % preallocation for unwanted conditions (not to be expanded)
on3 = cell(1,length(leftover));
for i = 1:length(leftover)                                                  % loop that makes a cell for all types of unwanted condition in the order you listed in input
    ind = cellfun(@(c) strcmp(c,leftover{i}),names);                                % that will be added to the end of the wanted cells
    inds2move = find(ind == 1); % indices where names match
    dur3{i} = durations{inds2move};
    on3{i} = onsets{inds2move};
end
        td= cat(2,dur2,dur3);
        to = cat(2,on2,on3);                                               % concatenation of expanded conditions to unwanted, compact conditions
        tn = cat(2,nam2,leftover);
        durations = td;
        onsets = to;        
        names = tn;
        save([strrep(f_path, '.mat', 'new.mat')], 'durations', 'names', 'onsets', '-mat');
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    