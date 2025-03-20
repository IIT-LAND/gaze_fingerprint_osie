function gaze_fingerprint_analysis_exp2_sess2(WEIGHT,PLOT,RUNONSERVER,NWORKERS,other_session,test)
% gaze_fingerprint_analysis
%
% Runs gaze fingerprint classifier analysis, optimized to run in parallel
%
%   INPUT
%
%   WEIGHT = 'none', 'duration','pupil'
%   PLOT = 1 for plotting fixation heatmaps, 0 if not
%   RUNONSERVER = 1 for running on server, else 0
%   NWORKERS = number of workers for parallel pool
%   other_session = session (other than 2) for which you want to compute
%   gaze classification
%   test = true if you want to test on nstim = 5, false if you want to run
%   on all 700 stim
%
%   Example usage:
%
%   WEIGHT = 'none'; PLOT = false; RUNONSERVER = false; NWORKERS = 5;
%   other_session = "session_02"; test = true;
%   gaze_fingerprint_analysis_suppl(WEIGHT, PLOT, RUNONSERVER, NWORKERS, other_session, test);
%

%% save input arguments to a structure
input_args = parse_input_args(WEIGHT, PLOT, RUNONSERVER, NWORKERS, other_session, test);

%% get experiment information
exp_info = get_exp_info(input_args);

%% path information
path_info = get_paths(input_args);

%% get pheno data
pheno_fname = fullfile(path_info.datapath,'pheno','pheno_OSIE_exp2.csv');
pheno_data = readtable(pheno_fname);

%% get information about target and comparison subjects
[target] = get_target_info(pheno_data);

%% get target fixation data
target.fixations = get_fix_data(target, pheno_data, path_info);

%% run gaze similarity analysis
gsa_res = gsa(target, input_args, exp_info, path_info);

%% write out results
output_results(gsa_res, path_info, exp_info);

%% plot heatmaps
if input_args.PLOT
    plot_heatmaps(target, input_args, exp_info, path_info)
end % if input_args.PLOT

end % function gaze_fingerprint_analysis


%% function to parse input arguments
function input_args = parse_input_args(WEIGHT, PLOT, RUNONSERVER, NWORKERS, other_session, test)

input_args.WEIGHT = WEIGHT;
input_args.PLOT = PLOT;
input_args.RUNONSERVER = RUNONSERVER;
input_args.NWORKERS = NWORKERS;
input_args.other_session = other_session;
input_args.test = test;

end % function parse_input_args

%% function to get paths
function path_info = get_paths(input_args)

if input_args.RUNONSERVER
  path_info.rootpath = '/media/DATA/Scratch/eye_tracking_data/gaze_fingerprinting_exp2';
else
  path_info.rootpath = '/Users/scrockford/Library/CloudStorage/OneDrive-FondazioneIstitutoItalianoTecnologia/gaze_fingerprinting_exp2';
end % if RUNONSERVER

path_info.codepath = fullfile(path_info.rootpath,'code');
path_info.datapath = fullfile(path_info.rootpath,'data');
path_info.resultpath = fullfile(path_info.rootpath,sprintf('resultstest', input_args.other_session));
path_info.plotpath = fullfile(path_info.rootpath,'plots');
path_info.stimpath = fullfile(path_info.rootpath,'stimuli');

end % function get_paths


%% function to grab experiment info
function exp_info = get_exp_info(input_args)

% set file stem
exp_info.fstem = sprintf('weight_%s',input_args.WEIGHT);

% set file stem
exp_info.plot_fstem = '.jpeg'; % '.pdf';

% smoothing radius for fixation heatmap
exp_info.smoothing_radius = 22;

% k nearest neighbors parameter
exp_info.k = 1;

if input_args.test 
    exp_info.nstim = 5; % test number of stimuli
else
    exp_info.nstim = 700; % true number of stimuli 
end % if test

exp_info.stim_dim = [600,800]; % stimulus dimensions

exp_info.entropy_nrow = 3; % addition for new way of calculating stationary entropy
exp_info.entropy_ncol = 4; % addition for new way of calculating stationary entropy

end % function get_exp_info

%% function to get info about target and comparison subjects
function [target] = get_target_info(pheno_data)

target.name = 'target';
target.mask = ismember(pheno_data.target_use,'Yes');
target.nsubs = sum(target.mask);
target.subids = pheno_data.subj_ID(target.mask);

end % function get_target_info

%% function to get fixation data
function fix_data = get_fix_data(target, pheno_data, path_info)

disp('Loading fixation data');

% pre-allocate empty dataframes
nsubs = target.nsubs;

fix_data.sess01{1, nsubs} = [];
fix_data.sess02{1, nsubs} = [];

for isub = 1:length(target.subids)

    sub_mask = ismember(pheno_data.subj_ID, target.subids{isub});
    sub_pheno_data = pheno_data(sub_mask,:);

    fpath2use = fullfile(path_info.datapath, ...
        sub_pheno_data.data_directory_name{1});

    % target session 1
    
    fpath2use_session1 = fullfile(fpath2use,sub_pheno_data.subj_ID{1}, ...
                                'session_01','fixations');
    fdir2use_session1 = dir(fullfile(fpath2use_session1,'*fixations.mat'));

    for file = fdir2use_session1'

        fname2use = fullfile(fpath2use_session1,file.name);
        disp(fname2use);

        if contains(fname2use, sub_pheno_data.subj_ID{1})
            tmp = load(fname2use);
            fix_data.sess01{isub} = tmp.data.fixations;
        else
            disp("no fixation data loaded")
        end % contains(fname2use, fnames.video_fname)
     end % file = fname2use_session1'

    % target session 2
    fpath2use_session2 = fullfile(fpath2use,sub_pheno_data.subj_ID{1}, ...
                                'session_02','fixations');
    fdir2use_session2 = dir(fullfile(fpath2use_session2,'*fixations.mat'));

    for file = fdir2use_session2'

        fname2use = fullfile(fpath2use_session2,file.name);
        disp(fname2use);

        if contains(fname2use, sub_pheno_data.subj_ID{1})
            tmp = load(fname2use);
            fix_data.sess02{isub} = tmp.data.fixations;
        else
            disp("no fixation data loaded")
        end % contains(fname2use, fnames.video_fname)

    end % file = fdir2use_session2'

end % for isub

end % function get_fix_data
%% gsa - main gaze similarity analysis function
function gsa_res = gsa(target, input_args, exp_info, path_info)

% pre-allocate arrays in memory
entropy_per_stim = zeros(target.nsubs, 2,exp_info.nstim);
intrasubject_r = zeros(target.nsubs, exp_info.nstim);
mean_intersubject_r = zeros(target.nsubs, exp_info.nstim);
sub_BY_stim_IDarray = zeros(target.nsubs, exp_info.nstim);
id_mat = nan(target.nsubs, target.nsubs, exp_info.nstim);
fingerprint_mat = nan(target.nsubs, target.nsubs, exp_info.nstim);
fingerprint_ratios = nan(target.nsubs, exp_info.nstim);
nfp_subs = nan(exp_info.nstim,1);
acc = nan(exp_info.nstim,1);
mean_fp_ratio = nan(exp_info.nstim,1);

tic;
parpool(input_args.NWORKERS);

parfor istim = 1:exp_info.nstim
%for istim = 1:nstim
    disp(fprintf('Stimulus %03d', istim));

    % compute target subs heatmaps ----------------------------------------
    hmap_res = get_heatmaps(target, istim, exp_info, input_args, path_info);

    % calculate intrasubject and intersubject similarity ------------------
    gsa_res_tmp = compute_gaze_similarity(hmap_res, target);

    % do kNN classifier and fingerprint ratio -----------------------------
    fp_res = fingerprint_classification(gsa_res_tmp);

    % put results together in these arrays --------------------------------
    id_mat(:,:,istim) = gsa_res_tmp.id_mat;
    fingerprint_mat(:,:,istim) = fp_res.fingerprint_mat;

    % grab correct or incorrect classification for each target subject
    sub_BY_stim_IDarray(:,istim) = fp_res.fingerprint_mat(logical(eye(size(fp_res.fingerprint_mat))));

    % grab entropy data
    entropy_per_stim(:,:,istim) = gsa_res_tmp.entropy;

    % grab intrasubject correlations
    intrasubject_r(:,istim) = gsa_res_tmp.intrasubject_r;

    % grab mean intersubject correlations
    mean_intersubject_r(:,istim) = gsa_res_tmp.mean_intersubject_r;

    % grab fingerprint ratios
    fingerprint_ratios(:,istim) = fp_res.fingerprint_ratios;

    % find number of subjects that can be fingerprinted
    nfp_subs(istim,1) = sum(sub_BY_stim_IDarray(:,istim), 'omitnan');

    % compute percent accuracy for fingerprinting using maximal correlation
    acc(istim,1) = nfp_subs(istim,1)./length(target.subids);

    % compute mean fingerprint_ratio across subjects
    mean_fp_ratio(istim,1) = mean(fp_res.fingerprint_ratios, 'omitnan');

end % for istim

% compute stationary entropy on fixation data -------------------------------
entropy_res = compute_stationaryentropy(target, exp_info);

% Pack everything into gsa_res
gsa_res.nstim = exp_info.nstim;
gsa_res.subids = target.subids;
gsa_res.mean_intrasubject_r = mean(intrasubject_r,1, 'omitnan');
gsa_res.intrasubject_r = intrasubject_r;
gsa_res.avg_mean_intersubject_r = mean(mean_intersubject_r,1, 'omitnan');
gsa_res.mean_intersubject_r = mean_intersubject_r;
gsa_res.sub_BY_stim_IDarray = sub_BY_stim_IDarray;
gsa_res.id_mat = id_mat;
gsa_res.fingerprint_mat = fingerprint_mat;
gsa_res.fingerprint_ratios = fingerprint_ratios;
gsa_res.fingerprintable_stimuli_per_sub = reshape(sum(sub_BY_stim_IDarray,2),target.nsubs,1);
gsa_res.nfp_subs = nfp_subs;
gsa_res.acc = acc;
gsa_res.mean_fp_ratio = mean_fp_ratio;
gsa_res.entropy_per_stim_sess1 = squeeze(entropy_per_stim(:,1,:));
gsa_res.entropy_per_stim_sess2 = squeeze(entropy_per_stim(:,2,:));
gsa_res.entropy = entropy_res.entropy;

toc;

end % function gsa


%% get_heatmaps function
function results = get_heatmaps(target, istim, exp_info, input_args, path_info)

nsubs = target.nsubs;
subids = target.subids;
stim_dim = exp_info.stim_dim;
smoothing_radius = exp_info.smoothing_radius;
WEIGHT = input_args.WEIGHT;
entropy_nrow = exp_info.entropy_nrow;
entropy_ncol = exp_info.entropy_ncol;

fixhmap1_tmp = zeros(nsubs,stim_dim(1)*stim_dim(2));
fixhmap2_tmp = zeros(nsubs,stim_dim(1)*stim_dim(2));
empty_hmap = zeros(nsubs,3);
res = zeros(stim_dim(1),stim_dim(2));

imgname = fullfile(path_info.stimpath,sprintf('1%03d.jpg',istim));

% compute target subs heatmaps ----------------------------------------
disp('Computing fixation heatmaps');
for isub = 1:nsubs
    % session 1
    fixationdata2use = target.fixations.sess01{isub}{istim};
    mask = isempty(fixationdata2use.fixX);
    if ~mask
        res = fixationHeatmap(fixationdata2use, imgname, smoothing_radius, 0, NaN, WEIGHT);
        fixhmap1_tmp(isub,:) = res(:);
    else
        fixhmap1_tmp(isub,:) = nan(1,size(fixhmap1_tmp,2));
        empty_hmap(isub,1) = 1;
    end % if ~mask

    % session 2
    if ~isempty(target.fixations.sess02{isub})
        fixationdata2use = target.fixations.sess02{isub}{istim};
        mask = isempty(fixationdata2use.fixX);
        if ~mask
            res = fixationHeatmap(fixationdata2use, imgname, smoothing_radius, 0, NaN, WEIGHT);
            fixhmap2_tmp(isub,:) = res(:);
        else
            fixhmap2_tmp(isub,:) = nan(1,size(fixhmap2_tmp,2));
            empty_hmap(isub,2) = 1;
        end % if ~mask
    end

end % for isub

results.istim = istim;
results.nsubs = nsubs;
results.subids = subids;
results.stim_img = imgname;
results.stim_dim = stim_dim;
results.fix_hmap1 = fixhmap1_tmp;
results.fix_hmap2 = fixhmap2_tmp;
results.empty_hmap = empty_hmap;
results.entropy_nrow = entropy_nrow;
results.entropy_ncol = entropy_ncol;

end % function get_heatmaps

%% compute_gaze_similarity function
% function results = compute_gaze_similarity(hmap_res)
function results = compute_gaze_similarity(hmap_res, target)

results.istim = hmap_res.istim;
results.nsubs = hmap_res.nsubs;
results.subids = hmap_res.subids;

results.id_mat = nan(results.nsubs,results.nsubs);
results.intrasubject_r = nan(results.nsubs,1);
results.intersubject_r = nan(results.nsubs,results.nsubs-1);
results.mean_intersubject_r = nan(results.nsubs,1);
results.entropy = nan(results.nsubs,2);

% calculate intrasubject correlation ----------------------------------
disp('Computing intrasubject correlations');
results.id_mat = corr(hmap_res.fix_hmap1',hmap_res.fix_hmap2');
mask = logical(eye(size(results.id_mat))); % extract the diagonal containg sub_i by sub_i values
results.intrasubject_r = results.id_mat(mask);

% % calculate entropy similarity ----------------------------------------
% disp('Computing entropy similarity');
% for isub = 1:length(results.subids)
%   tmp_hmap1 = reshape(hmap_res.fix_hmap1(isub,:), hmap_res.stim_dim(1), hmap_res.stim_dim(2));
%   tmp_hmap2 = reshape(hmap_res.fix_hmap2(isub,:), hmap_res.stim_dim(1), hmap_res.stim_dim(2));
%   results.entropy(isub,1) = entropy(tmp_hmap1);
%   results.entropy(isub,2) = entropy(tmp_hmap2);
% end

% *************************************************************************
% New addition for calculating entropy
disp('Computing entropy similarity');
total_rows = hmap_res.stim_dim(1);
total_cols = hmap_res.stim_dim(2);
entropy_nrow = hmap_res.entropy_nrow;
entropy_ncol = hmap_res.entropy_ncol;

for isub = 1:length(results.subids)
    tmp_fix1 = target.fixations.sess01{isub}{hmap_res.istim};
    tmp_fix2 = target.fixations.sess02{isub}{hmap_res.istim};
    results.entropy(isub,1) = stationary_entropy(tmp_fix1, entropy_nrow, entropy_ncol, total_rows, total_cols);
    results.entropy(isub,2) = stationary_entropy(tmp_fix2, entropy_nrow, entropy_ncol, total_rows, total_cols);
end % for isub
% *************************************************************************

% calculate intersubject correlation ----------------------------------
disp('Computing intersubject correlations');

for isub = 1:length(results.subids)
    
    % find subject with the maximal correlation with target
    tmp_vect = results.id_mat(isub,:);
    
    % compute fingerprint ratio from intra and intersubject correlations
    intersub_corrs = tmp_vect; intersub_corrs(isub) = []; % grab all intersubject correlations and pop-out the intrasubject correlation
    results.intersubject_r(isub,:) = intersub_corrs; % save intersubject correlations
    
    % save intersubject correlations
    tmp_mean_intersub_r = mean(intersub_corrs,'omitnan'); % compute mean intersubject correlation
    results.mean_intersubject_r(isub,:) = tmp_mean_intersub_r; % save mean intersubject correlation to results.intersubject_r
    
end % for isub


end % function compute_gaze_similarity


%% fingerprint_classification function
function fp_res = fingerprint_classification(gsa_res_tmp)

disp('Fingerprint classification');

fp_res.istim = gsa_res_tmp.istim;
fp_res.subids = gsa_res_tmp.subids;
fp_res.nsubs = gsa_res_tmp.nsubs;

% pre-allocate memory
fp_res.fingerprint_mat = zeros(fp_res.nsubs, fp_res.nsubs);
fp_res.fingerprint_ratios = nan(fp_res.nsubs, 1);

fp_res.final_r = [gsa_res_tmp.intrasubject_r, gsa_res_tmp.intersubject_r];

%===start new addition to fix bug 05.07.2023===============================

% the bug to fix is related to the fact that we have to take the id_mat in 
% gsa_res.id_mat and symmetrize it along the diagonal. Currently, that
% id_mat is not symmetrized along the diagonal. Instead, the upper triangle
% reflects session 1 target compared to session 2 distractors. The lower
% triangle has distractors from session 1. Because the intrasubject
% correlation is always session 1 of the target subject to session 2 of the
% target subject, we need to take the intersubject correlations computed
% when session 2 of the distractor subjects are used. Therefore, we need to
% take the id_mat and grab the upper triangle and replace the lower
% triangle with the values in the upper triangle. Once this is done, we can
% take that symmetrized id_mat and just grab each row and find the biggest.

% grab session 1 to session 2 id mat
id_mat = gsa_res_tmp.id_mat;

id_mat_symm = id_mat;
% tmp_mat = squeeze(id_mat(:,:,1));
mask = tril(id_mat)~=0;
tmp_mat = id_mat';
id_mat_symm(mask) = tmp_mat(mask); 
id_mat = id_mat_symm;

%===end new addition to fix bug 05.07.2023=================================

for isub = 1:fp_res.nsubs

    %===start new addition to fix bug 05.07.2023===========================
    
    % --- old code --------------------------------------------------------
    % find subject with the maximal correlation with target
    % tmp_vect = gsa_res_tmp.id_mat(isub,:);
    % --- old code --------------------------------------------------------

    % --- updated code fixing bug -----------------------------------------
    % find subject with the maximal correlation with target
    tmp_vect = id_mat(isub,:);
    % --- updated code fixing bug -----------------------------------------
    
    %===end new addition to fix bug 05.07.2023=============================

    % what is the max gaze similarity correlation
    max_corr = max(tmp_vect);

    % what is the target subject's intrasubject correlation
    curr_sub_intra = tmp_vect(isub);

    % is the target subject's intrasubject correlation the max correlation?
    fingerprint_result_curr_sub = curr_sub_intra==max_corr;
    
    % save that 0/1 value into fingerprint_mat
    fp_res.fingerprint_mat(isub,isub) = fingerprint_result_curr_sub; 

    % who was the winner subject (the one with the highest correlation)
    winner_sub = find(tmp_vect==max_corr);
    
    % fill in a 1 in fingerprint mat for the winner subject
    fp_res.fingerprint_mat(isub,winner_sub) = 1;
    fp_res.fingerprint_mat(winner_sub,isub) = 1;

    % compute fingerprint ratio from intra and intersubject correlations
    intersub_corrs = tmp_vect; intersub_corrs(isub) = []; % grab all intersubject correlations and pop-out the intrasubject correlation
    tmp_mean_intersub_r = mean(intersub_corrs,'omitnan'); % compute mean intersubject correlation
    fp_res.fingerprint_ratios(isub,1) = curr_sub_intra/tmp_mean_intersub_r; % compute fingerprint ratio and save into fingerprint_ratio field of fp_res

end % for isub

end % function fingerprint_classification


%% compute stationary entropy over all stimuli put together
function results = compute_stationaryentropy(target, exp_info)

nsubs = target.nsubs;
nstim = exp_info.nstim;
my_entropy = zeros(nsubs,2);

disp('Computing entropy on all stimuli put together');

for isub = 1:nsubs

    all_fix_sess01.fixX = [];
    all_fix_sess01.fixY = [];
    all_fix_sess02.fixX = [];
    all_fix_sess02.fixY = [];
    
    for istim = 1:nstim
        % concatenate all fixations across stimuli
        tmp_fix01 = target.fixations.sess01{isub}{istim};
        tmp_fix02 = target.fixations.sess02{isub}{istim};
    
        try
            all_fix_sess01.fixX = [all_fix_sess01.fixX, tmp_fix01.fixX];
            all_fix_sess01.fixY = [all_fix_sess01.fixY, tmp_fix01.fixY];
            all_fix_sess02.fixX = [all_fix_sess02.fixX, tmp_fix02.fixX];
            all_fix_sess02.fixY = [all_fix_sess02.fixY, tmp_fix02.fixY];
        catch
            all_fix_sess01.fixX = [all_fix_sess01.fixX, tmp_fix01.fix_x];
            all_fix_sess01.fixY = [all_fix_sess01.fixY, tmp_fix01.fix_y];
            all_fix_sess02.fixX = [all_fix_sess02.fixX, tmp_fix02.fix_x];
            all_fix_sess02.fixY = [all_fix_sess02.fixY, tmp_fix02.fix_y];
        end % try
    
    end % for istim

    total_rows = exp_info.stim_dim(1);
    total_cols = exp_info.stim_dim(2);
    entropy_nrow = exp_info.entropy_nrow;
    entropy_ncol = exp_info.entropy_ncol;
    my_entropy(isub, 1) = stationary_entropy(all_fix_sess01, entropy_nrow, entropy_ncol, total_rows, total_cols);
    my_entropy(isub, 2) = stationary_entropy(all_fix_sess02, entropy_nrow, entropy_ncol, total_rows, total_cols);

end % for isub

results.entropy = my_entropy;

end % function compute_stationaryentropy

%% stationary_entropy function
function result = stationary_entropy(fix_data, final_nrow, final_ncol, total_rows, total_cols)
% Compute stationary entropy on fixation data
%
% sge = stationary_entropy(fix_data, final_nrow, final_ncol, total_rows, total_cols)
%
% INPUT
%
%   fix_data = a structure with fixation data (e.g., fix_x, fix_y)
%   final_ncol = number of colummns you want in final state space 
%   final_nrow = number of rows you want in final state space
%   total_rows = number of rows in the actual stimulus (e.g., 600)
%   total_cols = number of columns in the actual stimulus (e.g., 800)
%
% Shannon's entropy: -1*sum_n(p*log2(p))
%
% Citation:
%
% Shiferaw, B., Downey, L., & Crewther, D. (2019). 
% A review of gaze entropy as a measure of visual scanning efficiency. 
% Neuroscience & Biobehavioral Reviews, 96, 353-366.


% calculate total number of fixations on the stimuli image
try
    nfix_tot = length(fix_data.fixX);
catch 
    nfix_tot = length(fix_data.fix_x);
end % try


% label state space grid cells
% stimulus_space = zeros(total_rows, total_cols);
% start_rows = 1:(total_rows/final_nrow):total_rows;
% start_columns = 1:(total_cols/final_ncol):total_cols;

% make a labeled state space rectangle for later counting up fixations
% within the labeled parts of state space
how_many_ss_rectangles = final_nrow * final_ncol;
how_many_rows_within_one_ss_rect = total_rows/final_nrow;
how_many_cols_within_one_ss_rect = total_cols/final_ncol;

my_counter = 0;
labeled_ss_rect = [];

for icol = 1:final_ncol
    tmp_column = [];

    for irow = 1:final_nrow
        my_counter = my_counter+1;
        tmp_mat = repmat(my_counter, how_many_rows_within_one_ss_rect, how_many_cols_within_one_ss_rect);
        tmp_column = [tmp_column; tmp_mat];
    end % for irow

    labeled_ss_rect = [labeled_ss_rect, tmp_column];

end % for icol


% find which cells in state space are the fixations
fix_count = zeros(how_many_ss_rectangles,1);

for ifix = 1:nfix_tot

    % grab current fixation coordinates
    try
        curr_fix_coords = [fix_data.fixY(ifix), fix_data.fixX(ifix)];
    catch
        curr_fix_coords = [fix_data.fix_y(ifix), fix_data.fix_x(ifix)];
    end

    % now figure out which cell in the state space grid the fixation falls
    % into
    idx2use = sub2ind(size(labeled_ss_rect), ...
        round(curr_fix_coords(1)), ...sure
        round(curr_fix_coords(2)));

    which_ss_cell = labeled_ss_rect(idx2use);
    fix_count(which_ss_cell) = fix_count(which_ss_cell)+1;
end % for ifix

% convert to proportion of fixations
prop_fix = fix_count./nfix_tot;

% remove zero entries in prop_fix 
prop_fix(prop_fix==0) = [];

%  calculate entropy
result = -1*sum(prop_fix .* log2(prop_fix));

end % function stationary_entropy



%% output_results function
function output_results(gsa_res, path_info, exp_info)

% write gsa_res to .mat file
save(fullfile(path_info.resultpath, sprintf('%s_gsa_res.mat',exp_info.fstem)), 'gsa_res');

% make up stimulus names for exporting files with column names per each
% stimulus
stim_names_str = cell(1, exp_info.nstim);
for istim = 1:gsa_res.nstim
    stim_names_str{istim} = sprintf('1%03d',istim);
end % for istim
stim_names = cell2table(stim_names_str','VariableNames',{'stimulus'});

% -------------------------------------------------------------------------
% table showing stats per stimulus
% save results into fingerprint results array
fingerprint_results(:,1) = gsa_res.nfp_subs;
fingerprint_results(:,2) = gsa_res.acc;
fingerprint_results(:,3) = gsa_res.mean_fp_ratio;
fingerprint_results(:,4) = gsa_res.mean_intrasubject_r;
fingerprint_results(:,5) = gsa_res.avg_mean_intersubject_r;

varNames = {'nsubs_fingerprint', 'acc_fingerprint', 'mean_fingerprint_ratio', ...
    'mean_intrasubject_correlation', 'avg_mean_intersubject_correlation'};
tab2write = [stim_names, array2table(fingerprint_results, 'VariableNames',varNames)];
fname2save = fullfile(path_info.resultpath,sprintf('%s_results_perStim.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% entropy session1 per subject x stimulus ------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.entropy_per_stim_sess1)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_entropy_session1_perSubStim.csv',exp_info.fstem));
writetable(tab2write, fname2save);

% entropy session1 per subject x stimulus ------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.entropy_per_stim_sess2)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_entropy_session2_perSubStim.csv',exp_info.fstem));
writetable(tab2write, fname2save);

% -------------------------------------------------------------------------
% intrasubject correlations per subject x stimulus ------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.intrasubject_r)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_intrasubjectR_perSubStim.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% -------------------------------------------------------------------------
% mean intersubject correlations per subject x stimulus -------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.mean_intersubject_r)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_mean_intersubjectR_perSubStim.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% -------------------------------------------------------------------------
% fingerprint ratios per subject x stimulus -------------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.fingerprint_ratios)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_fingerprintratios_perSubStim.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% -------------------------------------------------------------------------
% the classifier accuracy per subject x stimulus --------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.sub_BY_stim_IDarray)],'VariableNames',[{'subids'},stim_names_str]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_classification_accuracy_perSubStim.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% -------------------------------------------------------------------------
% the classification number of fingerprintable stimuli per subject --------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.fingerprintable_stimuli_per_sub)],'VariableNames',[{'subids','nfingerprintable_stimuli'}]);
fname2save = fullfile(path_info.resultpath, sprintf('%s_classification_nfingerprintableStimuli_perSub.csv',exp_info.fstem));
writetable(tab2write,fname2save);

% -------------------------------------------------------------------------
% entropy per subject and session -----------------------------------------
tab2write = cell2table([gsa_res.subids, num2cell(gsa_res.entropy)],'VariableNames',{'subids','entropy_session01','entropy_session02'});
fname2save = fullfile(path_info.resultpath, sprintf('%s_entropy_perSub.csv',exp_info.fstem));
writetable(tab2write,fname2save);

end % function output_results


%% plot_heatmaps function
function plot_heatmaps(target, input_args, exp_info, path_info)

% alpha level for transparency of heatmaps over stimuli
alpha_level = 0.75;

parfor istim = 1:exp_info.nstim

    % read in stimulus ------------------------------------------------
    imgname = fullfile(path_info.stimpath, sprintf('1%03d.jpg',istim));
    stim = imread(imgname);

    % make figure

    for isub = 1:length(target.subids)

        fig1 = figure; set(gcf,'color','w');

        % session 1
        fixationdata2use = target.fixations.sess01{isub}{istim};
        mask = isempty(fixationdata2use.fixX);
        if ~mask
            fixhmap1 = fixationHeatmap(fixationdata2use, imgname, exp_info.smoothing_radius, 0, NaN, input_args.WEIGHT);
        else
            fixhmap1 = zeros(exp_info.stim_dim(1),exp_info.stim_dim(2));
        end % if ~mask

        % session 2
        fixationdata2use = target.fixations.sess02{isub}{istim};
        mask = isempty(fixationdata2use.fixX);
        if ~mask
            fixhmap2 = fixationHeatmap(fixationdata2use, imgname, exp_info.smoothing_radius, 0, NaN, input_args.WEIGHT);
        else
            fixhmap2 = zeros(exp_info.stim_dim(1),exp_info.stim_dim(2));
        end % if ~mask

        subplot(1,2,1);
        imshow(stim); hold on; h = imshow(fixhmap1,'Colormap',parula);
        set(h, 'AlphaData', alpha_level);
        title(sprintf('%s session 1',target.subids{isub}));

        subplot(1,2,2);
        imshow(stim); hold on; h = imshow(fixhmap2,'Colormap',parula);
        set(h, 'AlphaData', alpha_level);
        title(sprintf('%s session 2',target.subids{isub}));

        path2use = fullfile(path_info.plotpath,'fixhmaps',input_args.WEIGHT, target.subids{isub});
        unix_str = sprintf('mkdir -p %s',path2use);
        unix(unix_str);

        fname2save = fullfile(path2use, sprintf('1%03d%s',istim,exp_info.plot_fstem));
        saveas(gcf, fname2save, exp_info.plot_fstem(2:end))
        close(fig1);

    end % for isub

end % for istim

end % function plot_heatmaps


%% fixationHeatmap function
function res = fixationHeatmap(data, imgname, smoothing_radius, PLOT, alpha_level,WEIGHT)
%
%   INPUT
%
%   data = fixation data - structure with
%               fix_x - set of x coordinates for each fixation
%               fix_y - set of y coordinates for each fixation
%               fix_duration - durations for each fixation
%
%   imgname = fullfile name to the stimulus
%
%   smoothing_radius = smoothing radius to use (e.g., 24)
%
%   PLOT = set to 1 to create the plot
%
%   alpha_level = alpha to set transparency. Set to NaN if not used.
%
%   WEIGHT = 'duration', 'pupil', or 'none'
%
%   OUTPUT
%
%   res = fixation heatmap
%

%% read in the stimulus
% get size of the stimulus which will constrain the size of the resulting heatmap
img = im2double(imread(imgname));
[h w ~] = size(img);
map = zeros([h w]);

% grab the fixation data
try
    data2useX = data.fixX;
    data2useY = data.fixY;
    dataDuration = data.fixDurations;
catch
    data2useX = data.fix_x;
    data2useY = data.fix_y;
    dataDuration = data.fix_duration;
end % try

% -------------------------------------------------------------------------
% remove out of bound coordinates
mask1 = data2useX>w | data2useX<1;
mask2 = data2useY>h | data2useY<1;
oob_mask = mask1 | mask2;
if sum(oob_mask)>0
    data2useX(oob_mask) = [];
    data2useY(oob_mask) = [];
    dataDuration(oob_mask) = [];
end % if sum(oob_mask)>0
% -------------------------------------------------------------------------

fix_x = max(1, min(round(data2useX), w));
fix_y = max(1, min(round(data2useY), h));
fix_duration = dataDuration;

% fill map with 1 where the fixations are
fix_x = floor(fix_x);
fix_y = floor(fix_y);
% loop over fixations
for k = 1:length(fix_x)
    if strcmp(WEIGHT,'duration')
        map(fix_y(k), fix_x(k)) = fix_duration(k);
    elseif strcmp(WEIGHT,'pupil')
        fix_pupil = data.fixPupilSize;
        map(fix_y(k), fix_x(k)) = fix_pupil(k);
    else
        map(fix_y(k), fix_x(k)) = 1;
    end % if strcmp(WEIGHT,'duration')
end % for k

% make the smoothing kernel
smoothingKernel = makeSmoothingKernel(smoothing_radius);

% smooth the fixation map
map = imfilter(map, smoothingKernel, 0);
% divide map by the max value in the map (this is what Xu et al., calls
% normalization
map = normalise(map);

res = map;

if PLOT
%     subplot(2,1,1);
%     imshow(img);
%     subplot(2,1,2);
    imshow(img); hold on;

    h = imshow(map,'Colormap',parula);
    if ~isnan(alpha_level)
        set(h,'AlphaData',alpha_level)
    end % if ~isnan(alpha_level)
end % if PLOT

end % function fixationHeatmap


%% normalise function
function [normalised] = normalise(map)
% [normalised] = normalise(map)
%

map = map - min(min(map));
s = max(max(map));
if s > 0
    normalised = map / s;
else
    normalised = map;
end

end % function normalise


%% function for making smoothing kernel
function res = makeSmoothingKernel(radius2use)
%
%   radius2use = 24
%

winSize = ceil(radius2use * 7);
res = fspecial('gaussian', [winSize winSize], radius2use);

end % function makeSmoothingKernel
