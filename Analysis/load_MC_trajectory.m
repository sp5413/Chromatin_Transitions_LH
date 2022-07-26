function [positions, LH_coordinates] = load_MC_trajectory(n_cores, n_linkers)

poc = 0;


 
links = n_linkers;

input_file_name = sprintf('trajectory_n800.txt', (n_cores)*n_linkers+n_cores);

input_file=fopen(input_file_name, 'r');
   
positions=fscanf(input_file,'%E %E %E\n',[3 inf]);

fclose(input_file);

   
LH_description_file = textread('LH_N0G6C22.in', '%s', 'delimiter', '\n');
LH_coordinates = zeros(28,3);

l = 1;
for i = [4:9 11:32]
   coords = strsplit(char(LH_description_file(i)));
   for j = 1:3
       LH_coordinates(l,j) = str2num(char(coords(j)));
   end
   l = l + 1;
end 
  
