function [all_lhs_for_trajectory, all_linkers_for_trajectory] = extract_coordinates(positions, cores, linkers, linker_histone, first, last, trajectory_step)

all_lhs_for_trajectory = double.empty;
all_linkers_for_trajectory = double.empty;

if(linker_histone)
   lines_per_frame = cores*4 + linkers*4 + cores*50 + cores*28     
else
   lines_per_frame = cores*4 + cores*50 + cores*28;
end


for frame_number = first:trajectory_step:last


    if mod(frame_number-1, 100)==0
        disp(sprintf('Frame number %d', frame_number))
    end

    frame = positions(:,(frame_number-1)*lines_per_frame+1:frame_number*lines_per_frame)';

    core_positions = zeros(cores,3);
    core_orientations = zeros(cores*3,3);
    core_index = 0;

    % ------------------- Histone cores ------------------------------------------------------
    %duz=length(type_file);
    for i = 1:cores        
        core_index = core_index + 1;
        core_positions(core_index,1:3) = frame((i-1)*4+1, 1:3);
        core_orientations((core_index-1)*3+1,1:3) = frame((i-1)*4+2, 1:3);
        core_orientations((core_index-1)*3+2,1:3) = frame((i-1)*4+3, 1:3);
        core_orientations((core_index-1)*3+3,1:3) = frame((i-1)*4+4, 1:3);        
    end
    % ---------------------------------------------------------------------------------------

    % ------------------ DNA linkers --------------------------------------------------------
    dna_linker_positions = zeros(linkers,3);
    dna_linker_orientations = zeros(linkers*3,3);
    dna_linker_index = 0;
    poz = cores*4;
    for i = 1:linkers
        dna_linker_index = dna_linker_index + 1;        
        dna_linker_positions(dna_linker_index,1:3) = frame(poz+(i-1)*4+1, 1:3);
        dna_linker_orientations((dna_linker_index-1)*3+1, 1:3) = frame(poz+(i-1)*4+2, 1:3);
        dna_linker_orientations((dna_linker_index-1)*3+2, 1:3) = frame(poz+(i-1)*4+3, 1:3);
        dna_linker_orientations((dna_linker_index-1)*3+3, 1:3) = frame(poz+(i-1)*4+4, 1:3);        
    end
    % ---------------------------------------------------------------------------------------
    
    links = linkers/(cores);
    dna = zeros(cores*2, links, 3);
    for i=1:cores
        
        fl_i = i - 1;
        sl_i = i;
        
        
        
        if(fl_i>0)
            index_in_dna = (i-1)*2+1;
            for j = 1:links
                hlp = dna_linker_positions((fl_i - 1)*links + j, :) - core_positions(i, :);
                hlp = hlp*inv(core_orientations((i-1)*3+1:(i-1)*3+3,1:3)); 
                dna(index_in_dna, j, :) = hlp';
               end
        end
        
        
        if(sl_i>0)
            index_in_dna = i*2;
            for j = 1:links
                hlp = dna_linker_positions((sl_i - 1)*links + j, :) - core_positions(i, :);
                hlp = hlp*inv(core_orientations((i-1)*3+1:(i-1)*3+3,1:3));
                dna(index_in_dna, j, :) = hlp';
                %index_in_dna = index_in_dna + 1;
            end
        end
                         
        
           
    end
    all_linkers_for_trajectory = [all_linkers_for_trajectory dna];

    % ----------------------------- Histone tails -----------------------  
    poz = cores*4;
    i   = linkers;   
    current_position = poz + (i-1)*4 + 4 + 1;   
    end_tail = (current_position - 1) + cores*50;      
    tails = zeros(cores*50,3);
    tail_counter = 0;
    for i = current_position:end_tail
        tail_counter = tail_counter + 1;
        tails(tail_counter,1) = frame(i,1);
        tails(tail_counter,2) = frame(i,2);
        tails(tail_counter,3) = frame(i,3);       
    end
    % ---------------------------------------------------------------------------------------

    if(linker_histone)       
        current_position = i;
        lh = zeros(cores, 28, 3);
        j = 0;       

        for i = 1:cores          
            lh(i,:,:) = frame(current_position + 1:current_position + 28,:);
            current_position = current_position + 28;
        end
    end

    for i = 1:cores
        for j = 1:28
            hlp = reshape(lh(i, j, :), 1, 3)  - core_positions(i, :);
            hlp = reshape(hlp,1,3)*inv(core_orientations((i-1)*3+1:(i-1)*3+3,1:3));
            lh(i, j, :) = hlp;
        end
    end

    all_lhs_for_trajectory = [all_lhs_for_trajectory lh];
end
