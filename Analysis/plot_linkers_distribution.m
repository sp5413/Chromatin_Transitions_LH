function [fig] = plot_linkers_distribution(cores, linkers, all_linkers, all_lhs, LH_coordinates, which_core, LH_length, fig, LH_1, LH_2, LH_3)


start_x = mean(LH_coordinates(1:3, 1));
start_y = mean(LH_coordinates(1:3, 2));
end_x = mean(LH_coordinates(end-2:end, 1));
end_y = mean(LH_coordinates(end-2:end, 2));

beads = linkers/(cores);

duz_whole = sqrt((start_x-end_x)^2 + (start_y-end_y)^2);
start_distance = sqrt((start_x)^2 + (start_y)^2);

angle = acos((abs(start_x-end_x)/duz_whole));




t = 0:0.01:(2*pi+0.02);
kx = 5.5*sin(t);
ky = 5.5*cos(t);


fl_x = -4.8:0.1:4.8;
fl_y = 2.5*ones(size(fl_x));

sl_x = -4.8:0.1:4.8;
sl_y = -2.5*ones(size(fl_x));

tl_y = -2.5:0.1:2.5;
tl_x = -4.8*ones(size(tl_y));

frl_y = -2.5:0.1:2.5;
frl_x =  4.8*ones(size(frl_y));

% hb1_x =  6.2; hb1_y = 0;
% hb2_x =  8.8; hb2_y = 0;
% hb3_x = 11.4; hb3_y = 0;


temp_all_linkers(:, :, 1) = cos(angle)*all_linkers(:, :, 1) - sin(angle)*all_linkers(:, :, 2);
temp_all_linkers(:, :, 2) = sin(angle)*all_linkers(:, :, 1) + cos(angle)*all_linkers(:, :, 2);
all_linkers(:, :, 1) = temp_all_linkers(:, :, 1);
all_linkers(:, :, 2) = temp_all_linkers(:, :, 2);


temp_all_lhs(:, :, 1) = cos(angle)*all_lhs(:, :, 1) - sin(angle)*all_lhs(:, :, 2);
temp_all_lhs(:, :, 2) = sin(angle)*all_lhs(:, :, 1) + cos(angle)*all_lhs(:, :, 2);
all_lhs(:, :, 1) = temp_all_lhs(:, :, 1);
all_lhs(:, :, 2) = temp_all_lhs(:, :, 2);


DNA_1 = int8.empty;
DNA_2 = int8.empty;
DNA_3 = int8.empty;

DNA_XY = zeros(50, 50);  % Histogram
DNA_XZ = zeros(50, 50);  % histogram

offset_1 = 14;
offset_2 = 25;
offset_3 = 25;


for i=1:(cores*2)
    for j=1:size(all_linkers,2)
        DNA_1 = [DNA_1 ceil(all_linkers(i, j, 1))];
        DNA_2 = [DNA_2 ceil(all_linkers(i, j, 2))];
        DNA_3 = [DNA_3 ceil(all_linkers(i, j, 3))];        
    end
end


for i=1:length(DNA_1)
    DNA_XY(DNA_1(i) + offset_1, DNA_2(i) + offset_2) = DNA_XY(DNA_1(i) + offset_1, DNA_2(i) + offset_2) + 1;
    DNA_XZ(DNA_1(i) + offset_1, DNA_3(i) + offset_3) = DNA_XZ(DNA_1(i) + offset_1, DNA_3(i) + offset_3) + 1;
end


DNA_XY = rot90(DNA_XY);
DNA_XY = flipud(DNA_XY);




linker_1 = zeros(cores, 2, beads);
linker_2 = zeros(cores, 2, beads);
linker_3 = zeros(cores, 2, beads);

which_cores = [which_core];

for i=which_cores
    l = 1;
    %for m=(i-1)*2+1:(i-1)*2+2 
    for m=(i-1)*2+1:i*2 
        m
        for j=1:beads
            linker_1(i, l, j) = mean(all_linkers(m, j:beads:end, 1));
            linker_2(i, l, j) = mean(all_linkers(m, j:beads:end, 2));
            linker_3(i, l, j) = mean(all_linkers(m, j:beads:end, 3));
        end
        l = l + 1;
    end
end

DNA_beads_1 = zeros(2,beads);
DNA_beads_2 = zeros(2,beads);
DNA_beads_3 = zeros(2,beads);

for j=1:beads
    %for i=1:cores
    for i = 1:which_cores
        DNA_beads_1(1, j) = DNA_beads_1(1, j) + linker_1(i, 1, j);
        DNA_beads_1(2, j) = DNA_beads_1(2, j) + linker_1(i, 2, j);

        DNA_beads_2(1, j) = DNA_beads_2(1, j) + linker_2(i, 1, j);
        DNA_beads_2(2, j) = DNA_beads_2(2, j) + linker_2(i, 2, j);
        
        DNA_beads_3(1, j) = DNA_beads_3(1, j) + linker_3(i, 1, j);
        DNA_beads_3(2, j) = DNA_beads_3(2, j) + linker_3(i, 2, j);        
    end    
end


DNA_beads_1 = DNA_beads_1/cores;
DNA_beads_2 = DNA_beads_2/cores;
DNA_beads_3 = DNA_beads_3/cores;

 



fig = fig + 1;
figp = figure(fig);    
subplot(2,1,1)
title(sprintf('Core %d',which_core))
hold on
%plot(all_linkers(:, 1:1:end, 1), all_linkers(:, 1:1:end, 2), 'r.')
for wi = which_core
    plot(all_linkers((wi*2-1):wi*2, 1:1:end, 1), all_linkers((wi*2-1):wi*2, 1:1:end, 2), 'r.')
end
    
hold on
plot(kx,ky,'k','LineWidth',2)

for j=1:6
    plot(LH_1(j),   LH_2(j),  'o', ...
        'MarkerEdgeColor',[210/255 105/255 30/255],...
        'MarkerFaceColor',[1.0 1.0 0.0],...
        'MarkerSize', 8)
end

for j=7:LH_length
    plot(LH_1(j), LH_2(j),  'o', ...
         'MarkerEdgeColor',[72/255  61/255 139/255],...
         'MarkerFaceColor',[00/255 206/255 209/255],...
         'MarkerSize', 7);
end

for wi = which_core
    plot(reshape(linker_1(wi, 1, :), 1, beads), reshape(linker_2(wi, 1, :), 1, beads), '.-b', 'LineWidth', 2);
    plot(reshape(linker_1(wi, 2, :), 1, beads), reshape(linker_2(wi, 2, :), 1, beads), '.-b', 'LineWidth', 2);
end

axis equal
axis([-10 30 -15 15])
set(gca,'fontsize',12)
LLX = xlabel('X coordinate');
set(LLX,'FontName','Arial');
set(LLX,'FontSize',12);        
LLY = ylabel('Y coordinate');
set(LLY,'FontName','Arial');
set(LLY,'FontSize',12); 
set(gca,'xtickmode', 'manual')
set(gca,'xtick',[-10:10:40])
set(gca,'xticklabel', {'-10', '0', '10', '20', '30', '40', '50'})
set(gca,'ytickmode', 'manual')
set(gca,'ytick',[-10:10:10])
set(gca,'yticklabel', {'-10', '0', '10', '15'})
plot(0,0,'o', ...
        'MarkerEdgeColor',[1.0 1.0 1.0],...
        'MarkerFaceColor',[1.0 1.0 1.0],...
        'MarkerSize', 11)


subplot(2,1,2) 
plot(all_linkers((which_core*2-1):which_core*2, 1:1:end, 1), all_linkers((which_core*2-1):which_core*2, 1:1:end, 3), 'r.')

hold on
plot(0,0,'o', ...
        'MarkerEdgeColor',[1.0 1.0 1.0],...
        'MarkerFaceColor',[1.0 1.0 1.0],...
        'MarkerSize', 11)
plot(fl_x,fl_y,'k','LineWidth',2)
plot(sl_x,sl_y,'k','LineWidth',2)
plot(tl_x,tl_y,'k','LineWidth',2)
plot(frl_x,frl_y,'k','LineWidth',2)



for wi = which_core
    plot(reshape(linker_1(wi, 1, :), 1, beads), reshape(linker_3(wi, 1, :), 1, beads), '.-b', 'LineWidth', 2);
    plot(reshape(linker_1(wi, 2, :), 1, beads), reshape(linker_3(wi, 2, :), 1, beads), '.-b', 'LineWidth', 2);
    

end
axis equal
%axis([-10 25 -15 15])
axis([-10 30 -15 15])
set(gca,'fontsize',12)
LLX = xlabel('X coordinate');
set(LLX,'FontName','Arial');
set(LLX,'FontSize',12);        
LLY = ylabel('Z coordinate');
set(LLY,'FontName','Arial');
set(LLY,'FontSize',12); 
set(gca,'xtickmode', 'manual')
set(gca,'xtick',[[-10:10:40]])
set(gca,'xticklabel', {'-10', '0', '10', '20', '30', '40', '50'})

set(gca,'ytickmode', 'manual')
set(gca,'ytick',[-10:10:10])
set(gca,'yticklabel', {'-10', '0', '10', '15'})
plot(0,0,'o', ...
        'MarkerEdgeColor',[1.0 1.0 1.0],...
        'MarkerFaceColor',[1.0 1.0 1.0],...
        'MarkerSize', 11)

r = 150;
new_name = sprintf('DNA_linker_fan_cores_%d_linkers_%d_core_%d.png',cores, linkers/cores, which_core);
print(figp,'-dpng',sprintf('-r%d', r), new_name);

return

        
end
