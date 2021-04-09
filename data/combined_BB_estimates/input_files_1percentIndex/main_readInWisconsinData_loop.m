function void = main_readInWisconsinData_loop(void)

clear all; close all; clc;

for th = 0.01:0.01:0.20
    var_calling_threshold = th;

    for i = 1:44
    
        [variants, n_variants, tp_weight] = GetIndexContactFreqs(i, var_calling_threshold);
        data(i).donor_iSNVs = variants(:,1);
        data(i).recipient_iSNVs = variants(:,2);
        data(i).recipient_reads = variants(:,3);
        data(i).recipient_reads_var = variants(:,4);
        data(i).n_variants_at_cutoff = n_variants;
        data(i).weight = tp_weight; % weight is 1 when contact and index are established; if bidirectional is considered, then weight is 1/2.
        
        clear vals

    end
    outfile = strcat('Wisconsin_inputData_', int2str(th*100), 'percent');
    
    save(outfile, 'data', 'var_calling_threshold');
    
    clear data
end
