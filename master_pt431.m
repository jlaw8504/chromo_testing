function master_pt431(seed, cond_num, steps, step_path, chromo_cmd)
cd(step_path);
%add condensin to the pT431_equilibrium.out file
add_condensin_ver3(seed, 'pT431_equilibrium.out', sprintf('pt431_s%d_c%d_0000.out', seed, cond_num), cond_num);
%define basename
basename = sprintf('pt431_s%d_c%d', seed, cond_num);
%Run first iteration
system(sprintf('%s -save %s_0000.out 5000 10 -continue',...
    chromo_cmd, basename))
%Save spring information
system(sprintf('grep spring %s_0000.out > springs_%s.txt',...
    basename, basename));
%% Start While loop
n = 1;
while n < steps + 1
    i = n - 1;
    if n < 10
        str_pad = '000';
    elseif n < 100
        str_pad = '00';
    elseif n < 1000
        str_pad = '0';
    else
        str_pad = [];
    end
    if i < 10
        str_padi = '000';
    elseif i < 100
        str_padi = '00';
    elseif i < 1000
        str_padi = '0';
    else
        str_padi = [];
    end
    %setup the filenames
    filename_i = sprintf('%s_%s%d.out', basename, str_padi, i);
    filename_n = sprintf('%s_%s%d.out', basename, str_pad, n);
    %Step the condensins
    condensin_step_ver5(filename_i, filename_n);
    %Run chromoShake
    system(sprintf('%s -save %s 5000 10 -continue', chromo_cmd, filename_n))
    %search for NaNs
    [~,result] = system(sprintf('grep nan %s | wc -l', filename_n));
    result_list = strsplit(result);
    nan_num = str2double(result_list{end-1});
    disp(strcat('Number of NaNs:',num2str(nan_num)));
    if nan_num > 1
        %remove the offending simulation iteration
        system(sprintf('rm %s', filename_n));
        [~,rseed_string] = system(sprintf('grep seed %s', filename_i));
        rseed_cell = strsplit(rseed_string);
        rseed = str2double(rseed_cell{end-1});
        new_rseed = rseed + 42;
        system(sprintf('sed -i ''s/random_number_seed %d/random_number_seed %d/g'' %s',...
            rseed, new_rseed, filename_i))
    else
        n = n + 1;
		system(sprintf('grep spring %s >> springs_%s.txt',...
            filename_n, basename));
		system(sprintf('rm %s', filename_i));
    end
end
