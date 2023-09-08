% bootstrap_all_accuracy.m
rootpath = '/Users/mlombardo/Dropbox/data/land_eye_tracking_data/gaze_fingerprinting';
load(fullfile(rootpath,'results','weight_none_gsa_res.mat'));

id_mat = gsa_res.id_mat; 
id_mat_symm = id_mat; 
tmp_mat = squeeze(id_mat(:,:,1)); 
mask = tril(tmp_mat)~=0; 
for istim = 1:size(id_mat,3)
    tmp_mat = squeeze(id_mat(:,:,istim)); 
    tmp_mat_rev = tmp_mat'; 
    tmp_mat2use = tmp_mat; 
    tmp_mat2use(mask) = tmp_mat_rev(mask); 
    id_mat_symm(:,:,istim) = tmp_mat2use; 
end
id_mat = id_mat_symm; 
tmp_ct = nan(size(id_mat,1),size(id_mat,2));

for isub = 1:size(id_mat,1)
    tmp_mat = squeeze(id_mat(isub,:,:))'; 
    tmp_ct(isub,:) = nanmedian(tmp_mat,1); 
end
data2use = tmp_ct;

na_mask = sum(isnan(gsa_res.intrasubject_r),2)>(700*0.15);
data2use = data2use(~na_mask,~na_mask);
subids = gsa_res.subids(~na_mask);

nboot = 10000;
for iboot = 1:nboot
    rng(iboot);
    boot_samples = randsample(1:length(subids), length(subids), true);
    subids_boot = subids(boot_samples);
    boot_mat = data2use(boot_samples, :);

    for isub = 1:size(boot_mat,1)
        max_sub = find(boot_mat(isub,:)==max(boot_mat(isub,:),[],'omitnan'));
        hit_miss(isub,1) = boot_samples(isub)==max_sub;
    end % for isub
    accuracy(iboot, 1) = sum(hit_miss)/size(boot_mat,1);
%     disp(sum(hit_miss)/size(boot_mat,1));
end % for iboot

% figure; hist(accuracy,100);
ci95 = prctile(accuracy,[2.5,97.5]);
tab2write = cell2table(num2cell([mean(accuracy),ci95]),'VariableNames',{'mean','low95','hi95'});
fname2save = fullfile(rootpath,'results','bootstrap_accuracy_95CIs.csv');
writetable(tab2write, fname2save);
