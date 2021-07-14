function SC_err_analysis(results_dir, dataset)

    load(sprintf('%s/SC/%s/Curvature.mat', results_dir, dataset), 'S', 'dS');
    load(sprintf('%s/SC/%s_Meta.mat', results_dir, dataset));
    
    rng(434983);
    
    N_control = 1000;
    which_K = [2, 5, 10, 50, 100, 250];
    neg_idx = arrayfun(@(i) randperm(length(S), N_control), [1:length(S)]', 'un', false);
    neg_idx = vertcat(neg_idx{:});
    
    [rho_pos, p_pos] = corr(S, S(knn(:,which_K)), 'rows', 'complete');
    [rho_neg, p_neg] = corr(S, S(neg_idx), 'rows', 'complete');
        
    frac_in_eb_pos = frac_in_95pct_CI(S(knn(:,which_K)), S, dS);
    frac_in_eb_neg = frac_in_95pct_CI(S(neg_idx), S, dS);
    
    SS2 = load(sprintf('%s/SC/%s_SS2/Curvature.mat', results_dir, dataset), 'S', 'dS');
    SS4 = load(sprintf('%s/SC/%s_SS4/Curvature.mat', results_dir, dataset), 'S', 'dS');

    frac_in_SS_pos(1) = frac_in_95pct_CI(S(ss_id{1}), SS2.S, SS2.dS);
    frac_in_SS_pos(2) = frac_in_95pct_CI(S(ss_id{2}), SS4.S, SS4.dS);
    
    neg_idx_SS2 = arrayfun(@(i) randperm(length(ss_id{1}))', [1:N_control], 'un', false);
    neg_idx_SS4 = arrayfun(@(i) randperm(length(ss_id{2}))', [1:N_control], 'un', false);
    neg_idx_SS2 = horzcat(neg_idx_SS2{:});
    neg_idx_SS4 = horzcat(neg_idx_SS4{:});
    
    frac_in_SS_neg(1,:) = frac_in_95pct_CI(S(ss_id{1}(neg_idx_SS2)), SS2.S, SS2.dS);
    frac_in_SS_neg(2,:) = frac_in_95pct_CI(S(ss_id{2}(neg_idx_SS4)), SS4.S, SS4.dS);

    save(sprintf('%s/SC/%s_ErrorAnalysis.mat', results_dir, dataset), 'which_K', 'rho_pos', 'p_pos', 'rho_neg', 'p_neg', ...
        'frac_in_eb_pos', 'frac_in_eb_neg', 'frac_in_SS_pos', 'frac_in_SS_neg');
    
    if (isequal(dataset, 'PBMC'))
        
        SS2 = load(sprintf('%s/SC/%s_SSTranscript2/Curvature.mat', results_dir, dataset), 'S', 'dS');
        SS4 = load(sprintf('%s/SC/%s_SSTranscript4/Curvature.mat', results_dir, dataset), 'S', 'dS');

        frac_in_SS_transcript_pos(1) = frac_in_95pct_CI(S(ss_transcript_id{1}), SS2.S, SS2.dS);
        frac_in_SS_transcript_pos(2) = frac_in_95pct_CI(S(ss_transcript_id{2}), SS4.S, SS4.dS);

        neg_idx_SS2 = arrayfun(@(i) randperm(length(ss_transcript_id{1}), N_control), [1:length(ss_transcript_id{1})]', 'un', false);
        neg_idx_SS4 = arrayfun(@(i) randperm(length(ss_transcript_id{2}), N_control), [1:length(ss_transcript_id{2})]', 'un', false);
        neg_idx_SS2 = vertcat(neg_idx_SS2{:});
        neg_idx_SS4 = vertcat(neg_idx_SS4{:});
        neg_idx_SS2 = ss_transcript_id{1}(neg_idx_SS2);
        neg_idx_SS4 = ss_transcript_id{2}(neg_idx_SS4);

        frac_in_SS_transcript_neg(1,:) = frac_in_95pct_CI(S(neg_idx_SS2), SS2.S, SS2.dS);
        frac_in_SS_transcript_neg(2,:) = frac_in_95pct_CI(S(neg_idx_SS4), SS4.S, SS4.dS);
        
        save(sprintf('%s/SC/%s_ErrorAnalysis.mat', results_dir, dataset), 'frac_in_SS_transcript_pos', 'frac_in_SS_transcript_neg', '-append');
        
    end
    
end

function f = frac_in_95pct_CI(val, S, dS)
    f = nanmean( ( val >= S-2*dS ) & ( val <= S+2*dS ) , 1);
end
