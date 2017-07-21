function condensin_step_ver5(infile,outfile)

% this code lets the condensin molecules on your DNA walk/extrude loops

% assign parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_sep = 1e-008; % standard distances betweem masses
spring_rest = 1e-008; % spring distance at rest
spring_const = 0.226195; % standard spring constant
spring_damp = 10; % how much factor by which the gammma spring is weaker than the normal spring constant
time = 'test'; % this is used to find the last time step
colors = 'test'; % used to assign mass colors
condensin_mass_color = 2; % color of condensin beads
DNA_mass_color = [1 4]; % color of DNA beads
mass_condensin = []; % this is used to look for condensin beads
springs_condensin = []; % this is used to look for springs bound to condensin
thresh = 3; % disatnce threshold for the condensin beads
spring_weak = 1000; % how much weaker the new alpha spring is

% open up the file
fid_in = fopen(infile);

% assign tline so that the lines can be looped through
tline = fgetl(fid_in);

m = 1; % used to assign mass coords
n = 1; % used to assign springs

% loop through all the lines to find the masses and appropriate springs
while ischar(tline)
    % find the mass coordinates
    if size(strfind(tline,'mass '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        % save the initial coordinates in a cell array
        mass_coords_in(m,1) = str2double(b{5});
        mass_coords_in(m,2) = str2double(b{6});
        mass_coords_in(m,3) = str2double(b{7});
        % record the color of each of the beads
        if size(b,2) == 7
            % rows without an eighth entry correpond to a color of 1 (red)
            mass_coords_in(m,4) = 1;
        elseif size(b,2) == 8
            % rows with an eighth entry have their color logged appropriately
            mass_coords_in(m,4) = str2double(b{8});
        else
        end
        mass_last(1,1) = str2double(b{3}); % assigns the last mass value in the file
        mass_coords_in(m,5) = str2double(b{3}); % puts the mass number in the table
        % increase the counter by 1
        m = m+1;
    else
    end
    if size(strfind(tline,'spring '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        
        % save the desired springs into a cell array
        springs(n,1) = str2double(b{3});
        springs(n,2) = str2double(b{4});
        springs(n,3) = str2double(b{6});
        
        % increase the counter by 1
        n = n+1;
    else
    end
    if size(strfind(tline,'Time '),1) ~= 0
        % assign the time line to the variable time, this will always pick the last one
        if strcmp(time,'test')==1
            time_init = tline;
        else
        end
        
        time = tline;
    else
    end
    if size(strfind(tline,'MassColors'),1) ~= 0
        % assign the color line to the variable colors
        colors = tline;
    else
    end
    % loop to the next line
    tline = fgetl(fid_in);
end

% close the file
fclose('all');

% this next part logs the colors of all the beads in order
if size(strfind(colors,'MassColors'),1) ~= 0
    % open up the file (again)
    fid_in = fopen(infile);
    % assign tline so that the lines can be looped through (starts it over)
    tline = fgetl(fid_in);
    while ischar(tline)
        % looks for the beginning of the last time step
        if strcmp(colors,tline) == 1
            % loop to the next line
            tline = fgetl(fid_in);
            % this assigns the coordinates (which are listed backwards for some reason)
            for d = 0:mass_last(1,1)
                % split the line into components
                % log all of the coordinates into a cell matrix
                mass_colors(mass_last(1,1)+1-d,1) = str2double(tline);
                % loop to the next line
                tline = fgetl(fid_in);
            end
        else
            % loop to the next line if it isn't the right one
            tline = fgetl(fid_in);
        end
    end
else
end

% close the file
fclose('all');

% if this is not the original file, it looks for the last time step
if size(strfind(time,'Time '),1) ~= 0
    % open up the file (again)
    fid_in = fopen(infile);
    % assign tline so that the lines can be looped through (starts it over)
    tline = fgetl(fid_in);
    % loop through all the lines
    while ischar(tline)
        % looks for the beginning of the last time step
        if strcmp(time,tline) == 1
            % loop to the next line
            tline = fgetl(fid_in);
            % this assigns the coordinates (which are listed backwards for some reason)
            for d = 0:mass_last(1,1)
                % split the line into components
                e = strsplit(tline);
                % log all of the coordinates into a cell matrix
                mass_coords(mass_last(1,1)+1-d,1) = str2double(e{1});
                mass_coords(mass_last(1,1)+1-d,2) = str2double(e{2});
                mass_coords(mass_last(1,1)+1-d,3) = str2double(e{3});
                mass_coords(mass_last(1,1)+1-d,4) = mass_colors(mass_last(1,1)+1-d,1);
                mass_coords(mass_last(1,1)+1-d,5) = mass_last(1,1)-d;
                % loop to the next line
                tline = fgetl(fid_in);
            end
        else
            % loop to the next line if it isn't the right one
            tline = fgetl(fid_in);
        end
    end
else
    % if this is the orginial document, then the coords were parsed earlier
    mass_coords = mass_coords_in;
end

% close the file
fclose('all');

% assign counters for the following section
m = 1; % used to assign DNA
n = 1; % used to assign condensin

% separate the condensin into a separate file
for z = 1:size(mass_coords)
    if max(mass_coords(z,4) == DNA_mass_color) == 1
        % assign everything with DNA colors to the DNA matrix
        mass_DNA(m,:) =  mass_coords(z,:);
        m = m+1;
    elseif max(mass_coords(z,4) == condensin_mass_color) == 1
        % assign everything with the condensin color to the condensin matrix
        mass_condensin(n,:) =  mass_coords(z,:);
        n = n+1;
    else
    end
end

m = 1; % used to assign condensin springs

% create a list of springs that contain condensin
if size(mass_condensin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,2) == mass_condensin(:,5)) == 1
            % if the spring contains a condensin bead, log it in a table
            springs_condensin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

m = 1; % reset counter to assign unbound DNA

% make a table of all unbound DNA

if size(springs_condensin,1)>0
    % if there is a condensin spring, select the DNA beads not bound to condensin
    for z = 1:size(mass_DNA,1)
        % loop through all the DNA
        if max(mass_DNA(z,5) == springs_condensin(:,1)) == 0
            % if this mass is not involved in any springs with condensin, log it
            mass_DNA_unb(m,:) = mass_DNA(z,:);
            m = m+1; % increase the counter by 1
        else
        end
    end
else
    % if there is no condensin, then all are unbound
    mass_DNA_unb = mass_DNA;
end

for z = 1:size(mass_condensin,1)/11
    % lets find all the A, B1, and B2 beads
    condensin_A(z,:) = mass_condensin(1+(11*(z-1)),:);
    condensin_B1(z,:) = mass_condensin(9+(11*(z-1)),:);
    condensin_B2(z,:) = mass_condensin(10+(11*(z-1)),:);
end

% preassign this with 0.5 so that it can be checked later
DNA_alpha = zeros(size(condensin_A,1),5)+0.5;
DNA_alpha_old = zeros(size(condensin_A,1),5)+0.5;

% use the mass numbers of the A, B1, and B2 beads to find out what they are bound to
for z = 1:size(condensin_A,1) % loop through the condensins
    for n = 1:size(springs_condensin,1) % loop through the springs
        % log all the DNA beads bound to A
        if condensin_A(z,5) == springs_condensin(n,2) &&...
                max(springs_condensin(n,1) == mass_condensin(:,5))==0 &&...
                springs_condensin(n,3) == spring_const;
            DNA_alpha(z,:) = mass_coords(springs_condensin(n,1)+1,:);
        else
        end
        if condensin_A(z,5) == springs_condensin(n,2) &&...
                max(springs_condensin(n,1) == mass_condensin(:,5))==0 &&...
                springs_condensin(n,3) < spring_const;
            DNA_alpha_old(z,:) = mass_coords(springs_condensin(n,1)+1,:);
        else
        end
        % log all the DNA beads bound to B1
        if condensin_B1(z,5) == springs_condensin(n,2) && max(springs_condensin(n,1) == mass_condensin(:,5))==0
            DNA_beta1(z,:) = mass_coords(springs_condensin(n,1)+1,:);
        else
        end
        % log all the DNA beads bound to B2
        if condensin_B2(z,5) == springs_condensin(n,2) && max(springs_condensin(n,1) == mass_condensin(:,5))==0
            DNA_beta2(z,:) = mass_coords(springs_condensin(n,1)+1,:);
        else
        end
    end
end

s = 1; % used to track the springs being deleted
y = 1; % used to track the springs being created
q = 2; % used to track the bound beads

% make empty matrices for a size check
spring_create = {};
spring_delete = {};

% assign a number that you know does not correspond to a bead for the check
bound_beads(1,1) = mass_last+1;

% now we're gonna loop through the condensin molecules and find out what lines need to be created/deleted
for z = 1:size(condensin_A,1)
    
    % calculate the distance between the beads adjacent to the two heads
    [sep_dist] = distance_between_3D_chromoshake(...
        [mass_coords(condensin_A(z,5)+2,1) mass_coords(condensin_B1(z,5),1)],...
        [mass_coords(condensin_A(z,5)+2,2) mass_coords(condensin_B1(z,5),2)],...
        [mass_coords(condensin_A(z,5)+2,3) mass_coords(condensin_B1(z,5),3)]);
    
    if DNA_alpha(z,5) ~= 0.5 && DNA_alpha_old(z,5) == 0.5 && sep_dist > thresh*mass_sep
        % case 1 corresponds to when the alpha spring is active and the condensin ends are far apart
        
        % this condition will make the alpha spring significantly weaker
        
        % log this spring to be deleted
        spring_delete{s,1} = sprintf('  spring %d %d %.1g %.6g',DNA_alpha(z,5),condensin_A(z,5),spring_rest,spring_const);
        s = s+1; % increase the counter by 1
        spring_create{y,1} = sprintf('  spring %d %d %.1g %.6g',DNA_alpha(z,5),condensin_A(z,5),spring_rest,spring_const/spring_weak);
        % increase the counters by 1
        y = y+1;
        q = q+1;
        
    elseif DNA_alpha(z,5) == 0.5 && DNA_alpha_old(z,5) ~= 0.5
        % case 2 corresponds to when there is only one weak alpha spring
        
        % this condition will attach the A end of condensin to the closest DNA bead
        
        m = 1; % counter for the beads that will be tested for distance
        clearvars alpha_bind_test % clear this variable
        alpha_bind_test = []; % sets the table for a size check
        
        % create a binding box for the DNA
        binding_box_x(z,1:2) = [condensin_A(z,1)-3*mass_sep condensin_A(z,1)+3*mass_sep];
        binding_box_y(z,1:2) = [condensin_A(z,2)-3*mass_sep condensin_A(z,2)+3*mass_sep];
        binding_box_z(z,1:2) = [condensin_A(z,3)-3*mass_sep condensin_A(z,3)+3*mass_sep];
        
        for n = 1:size(mass_DNA_unb) % loop through all the DNA beads
            if mass_DNA_unb(n,1) > binding_box_x(z,1) && mass_DNA_unb(n,1) < binding_box_x(z,2) &&...
                    mass_DNA_unb(n,2) > binding_box_y(z,1) && mass_DNA_unb(n,2) < binding_box_y(z,2) &&...
                    mass_DNA_unb(n,3) > binding_box_z(z,1) && mass_DNA_unb(n,3) < binding_box_z(z,2) &&...
                    max(mass_DNA_unb(n,5) == bound_beads(:))==0
                
                % log the close beads to be tested
                alpha_bind_test(m,:) = mass_DNA_unb(n,:);
                
                % increase the counter by 1
                m = m+1;
            else
            end
        end
        
        alpha_bind = [0.5 3*mass_sep]; % set a base parameter for comparison in the next for loop
        
        if size(alpha_bind_test,2)>0
            for n = 1:size(alpha_bind_test,1)
                % calculate the ditstance from A to this bead
                [sep_dist] = distance_between_3D_chromoshake(...
                    [alpha_bind_test(n,1) condensin_A(z,1)],...
                    [alpha_bind_test(n,2) condensin_A(z,2)],...
                    [alpha_bind_test(n,3) condensin_A(z,3)]);
                
                % reassign it if it is closer than the threshold
                if sep_dist <= alpha_bind(2)
                    alpha_bind = [alpha_bind_test(n,5) sep_dist];
                else
                end
            end
            
            if alpha_bind(1) ~= 0.5
                spring_delete{s,1} = sprintf('  spring %d %d %.1g %.6g',DNA_alpha_old(z,5),condensin_A(z,5),spring_rest,spring_const/spring_weak);
                s = s+1; % increase the counter by 1
                % create this new spring
                spring_create{y,1} = sprintf('  spring %d %d %.1g %.6g',alpha_bind(1),condensin_A(z,5),spring_rest,spring_const);
                % add this bead to the bound_beads table
                bound_beads(q,1) = alpha_bind(1);
                % increase the counters by 1
                y = y+1;
                q = q+1;
            else
            end
            
        else
        end
        
    elseif DNA_alpha(z,5) ~= 0.5 && DNA_alpha_old(z,5) == 0.5 && sep_dist <= thresh*mass_sep
        % case 3 corresponds to when the alpha spring is active and the condensin ends are close together
        
        % calculate the distance between the beads bound to A and B1
        [sep_dist_AB1] = distance_between_3D_chromoshake(...
            [DNA_alpha(z,1) DNA_beta1(z,1)],...
            [DNA_alpha(z,2) DNA_beta1(z,2)],...
            [DNA_alpha(z,3) DNA_beta1(z,3)]);
        
        % calculate the distance between the beads bound to A and B2
        [sep_dist_AB2] = distance_between_3D_chromoshake(...
            [DNA_alpha(z,1) DNA_beta2(z,1)],...
            [DNA_alpha(z,2) DNA_beta2(z,2)],...
            [DNA_alpha(z,3) DNA_beta2(z,3)]);
        
        % assign the B1 and B2 heads based on proximity to A
        if sep_dist_AB1 < sep_dist_AB2
            DNA_beta_close(z,:) = DNA_beta1(z,:);
            DNA_beta_far(z,:) = DNA_beta2(z,:);
            condensin_B_close(z,:) = condensin_B1(z,:);
        else
            DNA_beta_close(z,:) = DNA_beta2(z,:);
            DNA_beta_far(z,:) = DNA_beta1(z,:);
            condensin_B_close(z,:) = condensin_B2(z,:);
        end
        
        % find the vector between the two DNA beads that the B heads are bound to
        vect(z,1:3) = DNA_beta_far(z,1:3) - DNA_beta_close(z,1:3);
        
        % add the vector to the far bead to get a point around which to build the binding box
        test_spot(z,1:3) = DNA_beta_far(z,1:3) + vect(z,1:3);
        
        % create the region in which B will search for DNA to bind to
        binding_box_x(z,1:2) = [test_spot(z,1)-3*mass_sep test_spot(z,1)+3*mass_sep];
        binding_box_y(z,1:2) = [test_spot(z,2)-3*mass_sep test_spot(z,2)+3*mass_sep];
        binding_box_z(z,1:2) = [test_spot(z,3)-3*mass_sep test_spot(z,3)+3*mass_sep];
        
        m = 1; % counter for the beads that will be tested for distance
        clearvars beta_bind_test % clear this variable
        beta_bind_test = []; % sets the table for a size check
        
        % find the beads in the binding box
        for n = 1:size(mass_DNA_unb) % loop through all the DNA beads
            if mass_DNA_unb(n,1) > binding_box_x(z,1) && mass_DNA_unb(n,1) < binding_box_x(z,2) &&...
                    mass_DNA_unb(n,2) > binding_box_y(z,1) && mass_DNA_unb(n,2) < binding_box_y(z,2) &&...
                    mass_DNA_unb(n,3) > binding_box_z(z,1) && mass_DNA_unb(n,3) < binding_box_z(z,2) &&...
                    max(mass_DNA_unb(n,5) == bound_beads(:))==0
                
                % log the close beads to be tested
                beta_bind_test(m,:) = mass_DNA_unb(n,:);
                
                % increase the counter by 1
                m = m+1;
            else
            end
        end
        
        beta_bind = [0.5 2*mass_sep]; % set a base parameter for comparison in the next for loop
        
        if size(beta_bind_test,2)>0
            for n = 1:size(beta_bind_test,1)
                % calculate the ditstance from B to this bead
                [sep_dist] = distance_between_3D_chromoshake(...
                    [test_spot(z,1) beta_bind_test(n,1)],...
                    [test_spot(z,2) beta_bind_test(n,2)],...
                    [test_spot(z,3) beta_bind_test(n,3)]);
                
                % reassign it if it is closer than the threshold
                if sep_dist <= beta_bind(2)
                    beta_bind = [beta_bind_test(n,5) sep_dist];
                else
                end
            end
            
            % make sure a new distance was assigned
            if beta_bind(1) ~= 0.5
                % create this new spring
                spring_create{y,1} = sprintf('  spring %d %d %.1g %.6g',beta_bind(1),condensin_B_close(z,5),spring_rest,spring_const);
                % add this bead to the bound_beads table
                bound_beads(q,1) = beta_bind(1);
                % delete this spring
                spring_delete{s,1} = sprintf('  spring %d %d %.1g %.6g',DNA_beta_close(z,5),condensin_B_close(z,5),spring_rest,spring_const);
                % increase these counters by 1
                y = y+1;
                s = s+1;
                q = q+1;
            else
            end
            
        else
        end
        
    else
    end
end

% open up the file (again)
fid_in = fopen(infile);

% assign the name of the new file that is to be created
fid_out = fopen(outfile,'w');

% assign tline so that the lines can be looped through (starts it over)
tline = fgetl(fid_in);

% we're going to delete/create all the necessary springs
while ischar(tline)==1
    if size(strfind(tline,'spring '),1) ~= 0
        if size(spring_delete,1)>0 % if we're deleting springs
            
            m = 0; % set a base counter
            
            for z = 1:size(spring_delete,1)
                % check to see if the line matches any of the springs that will be deleted
                m = m + strcmp(tline,spring_delete{z,1});
            end
            
            if m == 0
                % reprint the line
                fprintf(fid_out,'%s\r\n',tline);
            else
            end
        else
            % reprint the line
            fprintf(fid_out,'%s\r\n',tline);
        end
    elseif strcmp(tline,'}')==1
        if size(spring_create)>0 %if we're creating springs
            for z = 1:size(spring_create)
                % loop through all the springs we're creating
                fprintf(fid_out,'%s\r\n',spring_create{z,1});
            end
            
            % ends the structure
            fprintf(fid_out,'}\r\n');
        else
            % otherwise, it reprints the line
            fprintf(fid_out,'%s\r\n',tline);
        end
    else
        % otherwise, it reprints the line
        fprintf(fid_out,'%s\r\n',tline);
    end
    
    % loop to the next line
    tline = fgetl(fid_in);
end

% close the files
fclose('all');

% time to itterate with ChromoShake!