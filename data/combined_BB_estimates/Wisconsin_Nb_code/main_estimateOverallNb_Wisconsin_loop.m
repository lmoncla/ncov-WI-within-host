function void = main_estimateOverallNb_Wisconsin_loop(void)

clear all; close all; clc;

th_cntr = 1;
for th = 0.01:0.01:0.20
    var_calling_threshold = th;

    infile = strcat('Wisconsin_inputData_', int2str(th*100), 'percent');
    title_str = strcat('Wisconsin cutoff: ', int2str(th*100), 'percent');
    outfile = strcat(infile, '_overallNb');

    load(infile);

    n_TPs = length(data);

    min_Nb = 1; max_Nb = 100;
    overall_Nb = min_Nb:max_Nb;
    
    figure(th_cntr); subplot(2,2,1); plot([0 1], [0 1], 'k--'); hold on; 
    plot([0 1], [var_calling_threshold var_calling_threshold], 'r--'); 
    plot([0 1], [(1-var_calling_threshold) (1-var_calling_threshold)], 'r--'); 
    xlabel('index SNV'); ylabel('contact SNV');
    
    subplot(2,2,2); xlabel('N_b'); ylabel('logL');
    
    for i = 1:n_TPs

        if data(i).n_variants_at_cutoff > 0
            [i data(i).n_variants_at_cutoff]
            data(i).logL_for_overall = GetLogL_forNb(data(i), var_calling_threshold, overall_Nb);
        else
            data(i).logL_for_overall = NaN*ones(1, length(min_Nb:max_Nb));
        end
        %figure(th_cntr); 
        %subplot(1,2,1); plot(data(i).donor_iSNVs, data(i).recipient_iSNVs, 'r.'); hold on; axis([0 1 0 1]);
        %subplot(1,2,2); plot(overall_Nb, data(i).logL_for_overall, 'r'); hold on;
    
        data(i).logL_for_overall(1:10)
        %pause
        
        %save(outfile, 'data', 'var_calling_threshold');
        subplot(2,2,1); plot(data(i).donor_iSNVs, data(i).recipient_iSNVs, 'b.'); hold on; 
        
        %xlabel('index SNV'); ylabel('contact SNV');

        subplot(2,2,2); plot(overall_Nb, data(i).logL_for_overall, 'b'); hold on;
        %xlabel('N_b'); ylabel('logL');
    end

    % zero-truncated Poisson distn
    cntr = 1;
    lambda_list = 0.1:0.1:50;
    for this_lambda = lambda_list
        pmf_poiss = poisspdf(overall_Nb, this_lambda)/(sum(poisspdf(overall_Nb, this_lambda)));
        for i = 1:n_TPs
            logL_allTPs(i) = GetLogLForMeanNb(pmf_poiss, overall_Nb, data(i).logL_for_overall);
            weights_allTPs(i) = data(i).weight;
            if isnan(logL_allTPs(i))
                logL_allTPs(i) = 0;
            end 
        end
    
        %overall_Nb_logL_allTPS(cntr) = sum(logL_allTPs);
        overall_Nb_logL_allTPs(cntr) = sum(logL_allTPs.*weights_allTPs); % weight TPs: if unidirectional, weight = 1; if bidirectional, weight = 1/2
        
        cntr = cntr + 1;
    end
    mean_Nb = (lambda_list)./(1-exp(-lambda_list));
    %figure(th_cntr); 
    subplot(2,2,3); plot(mean_Nb, overall_Nb_logL_allTPs);
    xlabel('mean N_b'); ylabel('logL');
    loc_max = find(overall_Nb_logL_allTPs == max(overall_Nb_logL_allTPs));
    
    locs_CI = find((max(overall_Nb_logL_allTPs) - overall_Nb_logL_allTPs) < 1.92)
    est_mean_Nb_allTPs = mean_Nb(loc_max)
    est_lambda_allTPs = lambda_list(loc_max)
    est_mean_Nb_allTPs_lowCI = mean_Nb(min(locs_CI));
    est_mean_Nb_allTPs_highCI = mean_Nb(max(locs_CI));
    
    est_pmf_poiss_allTPs = poisspdf(overall_Nb, est_lambda_allTPs)/(sum(poisspdf(overall_Nb, est_lambda_allTPs)));
    %figure(th_cntr); 
    subplot(2,2,4); bar(overall_Nb, est_pmf_poiss_allTPs)
    xlabel('transmission of N_b virions'); ylabel('probability');
 
    save(outfile, 'data', 'var_calling_threshold', 'overall_Nb', 'est_pmf_poiss_allTPs', 'est_mean_Nb_allTPs', 'est_mean_Nb_allTPs_lowCI', 'est_mean_Nb_allTPs_highCI');
    overall_Nb'
    est_pmf_poiss_allTPs'
    this_axis = axis; axis([0 50 this_axis(3) this_axis(4)]);
    
    subplot(2,2,1); title(title_str); 
    
    th_cntr = th_cntr + 1;
    
    clear data
    %pause
end
