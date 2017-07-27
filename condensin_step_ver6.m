function condensin_step_ver6(infile,outfile)
%% Isolate the initial mass coordinates and springs
[mass_coords, springs, hinges] = infile_mass_springs_id(infile);
%% Extract final coordinates from outfile
[mass_coords] = final_mass_coords(infile,mass_coords);
%% Main Step Function
stepfunction(mass_coords,springs,hinges,outfile);
end
