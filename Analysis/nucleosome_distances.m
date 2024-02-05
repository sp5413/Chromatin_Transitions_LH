function [frame_distances]=calculate_distances(core_positions)

duz = length(core_positions);

frame_distances = [];


for i = 6
    x1 = core_positions(i,1);
    y1 = core_positions(i,2);
    z1 = core_positions(i,3);
    
    x2 = core_positions(i+2,1);
    y2 = core_positions(i+2,2);
    z2 = core_positions(i+2,3);
    
    temp_distance = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    
    frame_distances = [frame_distances temp_distance];

    

end
