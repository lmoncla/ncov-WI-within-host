function [vals_pass, n_variants, tp_weight] = GetIndexContactFreqs(pairNumber, var_calling_threshold)

switch pairNumber
        case 1
            infile = 'pair_1a_beta-binomial-input.txt'; load(infile);
            vals = pair_1a_beta_binomial_input;
            tp_weight = 1/2;
        case 2
            infile = 'pair_1b_beta-binomial-input.txt'; load(infile);
            vals = pair_1b_beta_binomial_input;
            tp_weight = 1/2;
        case 3
            infile = 'pair_2a_beta-binomial-input.txt'; load(infile);
            vals = pair_2a_beta_binomial_input;
            tp_weight = 1/2;
        case 4
            infile = 'pair_2b_beta-binomial-input.txt'; load(infile);
            vals = pair_2b_beta_binomial_input;
            tp_weight = 1/2;
        case 5
            infile = 'pair_3a_beta-binomial-input.txt'; load(infile);
            vals = pair_3a_beta_binomial_input;
            tp_weight = 1/2;
        case 6
            infile = 'pair_3b_beta-binomial-input.txt'; load(infile);
            vals = pair_3b_beta_binomial_input;
            tp_weight = 1/2;
        case 7
            infile = 'pair_4_beta-binomial-input.txt'; load(infile);
            vals = pair_4_beta_binomial_input;
            tp_weight = 1;
        case 8
            infile = 'pair_5_beta-binomial-input.txt'; load(infile);
            vals = pair_5_beta_binomial_input;
            tp_weight = 1;
        case 9
            infile = 'pair_6a_beta-binomial-input.txt'; load(infile);
            vals = pair_6a_beta_binomial_input;
            tp_weight = 1/2;
        case 10
            infile = 'pair_6b_beta-binomial-input.txt'; load(infile);
            vals = pair_6b_beta_binomial_input;
            tp_weight = 1/2;
        case 11
            infile = 'pair_7_beta-binomial-input.txt'; load(infile);
            vals = pair_7_beta_binomial_input;
            tp_weight = 1;
        case 12
            infile = 'pair_8a_beta-binomial-input.txt'; load(infile);
            vals = pair_8a_beta_binomial_input;
            tp_weight = 1/2;
        case 13
            infile = 'pair_8b_beta-binomial-input.txt'; load(infile);
            vals = pair_8b_beta_binomial_input;
            tp_weight = 1/2;
        case 14
            infile = 'pair_9a_beta-binomial-input.txt'; load(infile);
            vals = pair_9a_beta_binomial_input;
            tp_weight = 1/2;
        case 15
            infile = 'pair_9b_beta-binomial-input.txt'; load(infile);
            vals = pair_9b_beta_binomial_input;
            tp_weight = 1/2;
        case 16
            infile = 'pair_10a_beta-binomial-input.txt'; load(infile);
            vals = pair_10a_beta_binomial_input;
            tp_weight = 1/2;
        case 17
            infile = 'pair_10b_beta-binomial-input.txt'; load(infile);
            vals = pair_10b_beta_binomial_input;
            tp_weight = 1/2;
        case 18
            infile = 'pair_11a_beta-binomial-input.txt'; load(infile);
            vals = pair_11a_beta_binomial_input;
            tp_weight = 1/2;
        case 19
            infile = 'pair_11b_beta-binomial-input.txt'; load(infile);
            vals = pair_11b_beta_binomial_input;
            tp_weight = 1/2;
        case 20
            infile = 'pair_12a_beta-binomial-input.txt'; load(infile);
            vals = pair_12a_beta_binomial_input;
            tp_weight = 1/2;
        case 21
            infile = 'pair_12b_beta-binomial-input.txt'; load(infile);
            vals = pair_12b_beta_binomial_input;
            tp_weight = 1/2;
        case 22
            infile = 'pair_13a_beta-binomial-input.txt'; load(infile);
            vals = pair_13a_beta_binomial_input;
            tp_weight = 1/2;
        case 23
            infile = 'pair_13b_beta-binomial-input.txt'; load(infile);
            vals = pair_13b_beta_binomial_input;
            tp_weight = 1/2;
        case 24
            infile = 'pair_14_beta-binomial-input.txt'; load(infile);
            vals = pair_14_beta_binomial_input;
            tp_weight = 1;
        case 25
            infile = 'pair_15_beta-binomial-input.txt'; load(infile);
            vals = pair_15_beta_binomial_input;
            tp_weight = 1;
        case 26
            infile = 'pair_16_beta-binomial-input.txt'; load(infile);
            vals = pair_16_beta_binomial_input;
            tp_weight = 1;
        case 27
            infile = 'pair_17_beta-binomial-input.txt'; load(infile);
            vals = pair_17_beta_binomial_input;
            tp_weight = 1;
        case 28
            infile = 'pair_18_beta-binomial-input.txt'; load(infile);
            vals = pair_18_beta_binomial_input;
            tp_weight = 1;
        case 29
            infile = 'pair_19_beta-binomial-input.txt'; load(infile);
            vals = pair_19_beta_binomial_input;
            tp_weight = 1;
        case 30
            infile = 'pair_20a_beta-binomial-input.txt'; load(infile);
            vals = pair_20a_beta_binomial_input;
            tp_weight = 1/2;
        case 31
            infile = 'pair_20b_beta-binomial-input.txt'; load(infile);
            vals = pair_20b_beta_binomial_input;
            tp_weight = 1/2;
        case 32
            infile = 'pair_21a_beta-binomial-input.txt'; load(infile);
            vals = pair_21a_beta_binomial_input;
            tp_weight = 1/2;
        case 33
            infile = 'pair_21b_beta-binomial-input.txt'; load(infile);
            vals = pair_21b_beta_binomial_input;
            tp_weight = 1/2;
        case 34
            infile = 'pair_22a_beta-binomial-input.txt'; load(infile);
            vals = pair_22a_beta_binomial_input;
            tp_weight = 1/2;
        case 35
            infile = 'pair_22b_beta-binomial-input.txt'; load(infile);
            vals = pair_22b_beta_binomial_input;
            tp_weight = 1/2;
        case 36
            infile = 'pair_23_beta-binomial-input.txt'; load(infile);
            vals = pair_23_beta_binomial_input;
            tp_weight = 1;
        case 37
            infile = 'pair_24_beta-binomial-input.txt'; load(infile);
            vals = pair_24_beta_binomial_input;
            tp_weight = 1;
        case 38
            infile = 'pair_25a_beta-binomial-input.txt'; load(infile);
            vals = pair_25a_beta_binomial_input;
            tp_weight = 1/2;
        case 39
            infile = 'pair_25b_beta-binomial-input.txt'; load(infile);
            vals = pair_25b_beta_binomial_input;
            tp_weight = 1/2;
        case 40
            infile = 'pair_26a_beta-binomial-input.txt'; load(infile);
            vals = pair_26a_beta_binomial_input;
            tp_weight = 1/2;
        case 41
            infile = 'pair_26b_beta-binomial-input.txt'; load(infile);
            vals = pair_26b_beta_binomial_input;
            tp_weight = 1/2;
        case 42
            infile = 'pair_27a_beta-binomial-input.txt'; load(infile);
            vals = pair_27a_beta_binomial_input;
            tp_weight = 1/2;
        case 43
            infile = 'pair_27b_beta-binomial-input.txt'; load(infile);
            vals = pair_27b_beta_binomial_input;
            tp_weight = 1/2;
        case 44
            infile = 'pair_28_beta-binomial-input.txt'; load(infile);
            vals = pair_28_beta_binomial_input;
            tp_weight = 1;
        otherwise
            error('not valid')
end

locs_pass = intersect(find(vals(:,1) > var_calling_threshold), find(vals(:,1) < (1 - var_calling_threshold)));
vals_pass = vals(locs_pass,:);

locs_contact_below = find(vals_pass(:,2) < var_calling_threshold);
vals_pass(locs_contact_below, 2) = 0;
locs_contact_above = find(vals_pass(:,2) > (1 - var_calling_threshold));
vals_pass(locs_contact_above, 2) = 1;

n_variants = length(locs_pass);
