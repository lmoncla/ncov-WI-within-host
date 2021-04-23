function void = main_PlotOverallEstimateNb_Wisconsin_loop(void)

clear all; close all; clc;

vals = [];
vals_low = [];
vals_high = [];

thlist = 0.01:0.01:0.20;
for th = thlist
    var_calling_threshold = th;

    infile = strcat('Wisconsin_inputData_', int2str(th*100), 'percent_overallNb');

    load(infile);

    n_TPs = length(data);

    %plot(th, est_mean_Nb_allTPS, 'ko'); hold on;
    
    vals = [vals est_mean_Nb_allTPs];
    vals_low = [vals_low est_mean_Nb_allTPs_lowCI];
    vals_high = [vals_high est_mean_Nb_allTPs_highCI];
    
end

% plot(th, est_mean_Nb_allTPS, 'ko'); hold on;
plot(thlist, vals, 'k'); hold on; plot(thlist, vals, 'ko'); hold on;
plot(thlist, vals_low, 'r'); hold on; plot(thlist, vals_low, 'ro'); hold on;
plot(thlist, vals_high, 'r'); hold on; plot(thlist, vals_high, 'ro'); hold on;
xlabel('variant calling threshold');
ylabel('max. likelihood estimate of mean Nb');  

plot(thlist, ones(size(thlist))*1.21, 'm--');
    