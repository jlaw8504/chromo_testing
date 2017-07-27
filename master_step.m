function master_step(seed, infile, basename, step_path, chromo_cmd, steps_per_output,output_num,max_steps,cond_num,is_this_continutation,current_step)
%% Initial setup
cd(step_path); %Set initial path
if is_this_continutation == 0 %If this starting anew
    master_outfile = sprintf('%s.out',basename); %Create an the master .out file name
    add_condensin_ver3(seed, infile, master_outfile, cond_num); %Add condensins
    system(sprintf('%s -save %s %d %d -continue',...
        chromo_cmd, master_outfile,steps_per_output,output_num)); %Run first set of iterations
else % If this is continuing from another .out file
    master_outfile = infile;
end
system(sprintf('grep spring %s > springs_%s.txt',...
    master_outfile, basename)); % Save the initial spring information
%% While loop
n = current_step+1; %Create a counter at the new step (1 if anew)
input = master_outfile; % Set the first input file to be the master_outfile
while n <= max_steps %Run simulation to maximum step number
    if n < 10 % Create a number padding for the string
        str_pad = '000';
    elseif n < 100
        str_pad = '00';
    elseif n < 1000
        str_pad = '0';
    else
        str_pad = [];
    end
    condensin_step_ver6(input,'temp.cfg'); %Create a new .cfg with new step information
    system(sprintf('%s -save temp.out %d %d temp.cfg',...
        chromo_cmd,steps_per_output,output_num)); %Run interations on temp.cfg
    [~,output] = system('grep -n -B2 -m1 "Time " temp.out | head -n1 | cut -d"-" -f1'); %Find the linenumber of the first timestep to isolate
    cut_off_point = str2double(output); %Convert the linenumber from string to double
    system(sprintf('sed "1,%dd" temp.out >> %s',cut_off_point,master_outfile)); %Remove the region before the beginning of the timesteps and added that to the end of master_outfile
    name_update = sprintf('%s_%s%d.out',basename,str_pad,n); %Create a new string that has the updated name (to easily track status)
    system(sprintf('rename %s %s', master_outfile, name_update)); % Rename the master .out to the new, updated name
    master_outfile = name_update; %Set the new updated name to the master_outfile variable so it can be properly referenced
    input = 'temp.out'; %Set temp.out to input so that it can re-enter the cycle of creating a new .cfg and so on
    n = n+1; %Increase the counter
    system(sprintf('grep spring temp.cfg >> springs_%s.txt', basename)); %Output the new spring information in the temp.cfg file to the main springs file
end
end