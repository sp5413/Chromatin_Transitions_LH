fig = 0;
close all force
clc


load_all = 1;


if(load_all)   
   clear all
   fig = 0;
   

   LH_length = 28;
   cores = 100;
   beads = 7;
   


   linker_histone = 1;
   linkers = beads*cores;
   linkers = beads*(cores);
   links = beads;
   
    
   
   all_lhs     = double.empty;
   all_linkers = double.empty;
   
   SETS = [1];

       
       fig=0;
              

               [positions, LH_coordinates] = load_MC_trajectory(cores, beads);
               
               if(linker_histone)
                   lines_per_frame = cores*4 + linkers*4 + cores*50 + cores*28     
               else
                   lines_per_frame = cores*4 + cores*50 + cores*28;
               end

               number_of_frames = floor(length(positions)/lines_per_frame);



               first = 1;
               last = number_of_frames;
               trajectory_step = 1;

                 [all_lhs_for_trajectory, all_linkers_for_trajectory] = extract_coordinates(positions, cores, linkers, linker_histone, first, last, trajectory_step);
                 all_lhs     = [all_lhs all_lhs_for_trajectory];
                 all_linkers = [all_linkers all_linkers_for_trajectory]


end  

size(all_lhs);

clear positions



which_core = 90;
disp('plot_linkers_distribution(cores, linkers, all_linkers, all_lhs, LH_coordinates, which_core, LH_length, fig, LH_1, LH_2, LH_3');
fig = plot_linkers_distribution(cores, linkers, all_linkers, all_lhs, LH_coordinates, which_core, LH_length, fig, LH_1, LH_2, LH_3);
disp('-------------------------------------------------------------------')
