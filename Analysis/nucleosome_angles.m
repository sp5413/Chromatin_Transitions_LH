function [triplet_angles]=calculate_triplet_angles_consecutive(core_positions)

duz = length(core_positions);

triplet_angles = [];

for i = 6

    x1 = core_positions(i,1);
    y1 = core_positions(i,2);
    z1 = core_positions(i,3);

    x2 = core_positions(i+1,1);
    y2 = core_positions(i+1,2);
    z2 = core_positions(i+1,3);

    x3 = core_positions(i+2,1);
    y3 = core_positions(i+2,2);
    z3 = core_positions(i+2,3);

    v12x = x1 - x2;
    v12y = y1 - y2;
    v12z = z1 - z2;

    v23x = x3 - x2;
    v23y = y3 - y2;
    v23z = z3 - z2;

    v32x = x2 - x3;
    v32y = y2 - y3;
    v32z = z2 - z3;


    V12 = [v12x v12y v12z];
    V23 = [v23x v23y v23z];
    V32 = [v32x v32y v32z];


    temp_angle = acos(dot(V12,V23)/(sqrt(dot(V12,V12))*sqrt(dot(V23,V23))))*180/pi;

          
    triplet_angles = [triplet_angles temp_angle];
end
