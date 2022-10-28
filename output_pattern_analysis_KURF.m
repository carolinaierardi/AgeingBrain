%Script Name: output_pattern_analysis.m
%Date: 16/10/2022
%Author: CMI
%Version: 3.0
%Purpose: output pattern similarity after learning for young and old,
%duplicate analysis
%Notes: we will look at the similarity in output patterns for young and old
%learning as well as analyse the input pattern similarity for duplicate
%patterns in each case

close all
clear all
clc

%% Output pattern young

young_nodes_fire_normal = readmatrix('Young_OP_normal.csv');
valid_input_patterns = readmatrix('IP_after_threshold_normal.csv');

nodes_fired_y = readmatrix('Young_nodes_fired.csv'); %get the output patterns for young learning
nodes_fired_sorted_y = sort(nodes_fired_y); %put them in order

similar6_y = zeros(1,length(nodes_fired_sorted_y(1,:))); %create a matrix to store how many patterns were duplicates

for ii = 1:length(nodes_fired_sorted_y(1,:)) %for all the output patterns

    %verify if the each node (in each position) firing is the same 
condition1 = find(nodes_fired_sorted_y(1,:) == nodes_fired_sorted_y(1,ii));
condition2 = find(nodes_fired_sorted_y(2,:) == nodes_fired_sorted_y(2,ii));
condition3 = find(nodes_fired_sorted_y(3,:) == nodes_fired_sorted_y(3,ii));
condition4 = find(nodes_fired_sorted_y(4,:) == nodes_fired_sorted_y(4,ii));
condition5 = find(nodes_fired_sorted_y(5,:) == nodes_fired_sorted_y(5,ii));
condition6 = find(nodes_fired_sorted_y(6,:) == nodes_fired_sorted_y(6,ii));
    
    %find whether they overlap across nodes
intersect1 = intersect(condition1,condition2);
intersect2 = intersect(intersect1,condition3);
intersect3 = intersect(intersect2,condition4);
intersect4 = intersect(intersect3,condition5);
intersect5 = intersect(intersect4,condition6);


    %if all the conditions are met, add the number to the matrix, -1 to
    %account for the patterns comparing itself
similar6_y(ii) = length(intersect5) - 1;

end

length(find(similar6_y)) %how many are duplicates
max(similar6_y) %what is the maximum number of duplicates


% now, we want to test for the similarity of input patterns in the output
% patterns that are duplicates

%initalise loop variable
r = 0;

duplicates_y = zeros; %set up matrix of zeros

for o = 1:length(nodes_fired_sorted_y) %compare all the input patterns
    for y = 1:length(nodes_fired_sorted_y) %with all other input patterns

        if o == y %don't compare if they're the same
            continue
        end
        
        equality = intersect(nodes_fired_sorted_y(:,o), nodes_fired_sorted_y(:,y)); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 6 %if they are duplicates or unique, continue to next one
                disp(o) %show the index of the old OP
                disp(y) %show the index of the young OP
                r = r + 1; %adjust indexing
                duplicates_y(r,1) = o; %first column shows first OP
                duplicates_y(r,2) = y; %second column shows second OP

        else %otherwise add to matrix in the row of the common nodes
        continue
        end
    end
    if length(duplicates_y(:,1)) == round(length(find(similar6_y))*0.01) %we only need a 1% sample
        break %so stop the loop when we find it
    end
end

%once I've found the duplicates, I will compare their input patterns

for d = 1:length(duplicates_y(:,1))
    first = valid_input_patterns(:,duplicates_y(d,1));
    second = valid_input_patterns(:,duplicates_y(d,2));
    inters = intersect(find(first), find(second));
    sim_dup_y(d) = length(inters);
end

%Uncomment the next line to save matrix as a file (and not have to run this
%again)
% writematrix(sim_dup_y, 'Young_identical_sim.csv')

%% unique patterns

unique_y = zeros(1, length(nodes_fired_sorted_y(1,:))); %set a matrix to zeros

for ee = 1:length(nodes_fired_sorted_y(1,:)) %for all the input patterns

    %check the nodes firing in the same position
    condition1 = find(nodes_fired_sorted_y(1,:) == nodes_fired_sorted_y(1,ee));
    condition2 = find(nodes_fired_sorted_y(2,:) == nodes_fired_sorted_y(2,ee));
    condition3 = find(nodes_fired_sorted_y(3,:) == nodes_fired_sorted_y(3,ee));
    condition4 = find(nodes_fired_sorted_y(4,:) == nodes_fired_sorted_y(4,ee));
    condition5 = find(nodes_fired_sorted_y(5,:) == nodes_fired_sorted_y(5,ee));
    condition6 = find(nodes_fired_sorted_y(6,:) == nodes_fired_sorted_y(6,ee));

    %if none of these conditions are met, add it to the matrix
    if length(condition1) == 0 && length(condition2) == 0 && length(condition3) == 0 && length(condition4) == 0 && length(condition5) == 0 && length(condition6) == 0
        unique_y(ee) = 1;
    end
end

find(unique_y) %how many unique patterns are there?

%% Using a representative sample, checking node similarity for 5-1 nodes in common

%now, to determine how many patterns have 5-1 node in common, we must pick
%a representative sample to compare the nodes to.

sample = round(length(nodes_fired_sorted_y(1,:))*0.01);
%we will select a sample of 10%
randomselection = randi(length(nodes_fired_sorted_y(1,:)), [1,sample]);
interest_inputs_y = nodes_fired_sorted_y(:,randomselection); %get the patterns for the randomly selected values 

similarity_values = [1:5]; %values I want to test for similarity
similarity_analysis_y = zeros(length(similarity_values),length(interest_inputs_y)); 
                                    %matrix with values of interest and
                                    %patterns of interest

for i = 1:length(interest_inputs_y) %random selection
    for ii = 1:length(nodes_fired_sorted_y) %all the patterns
        
        equality = intersect(interest_inputs_y(:,i), nodes_fired_sorted_y(:,ii)); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 0 || commonnodes == 6 %if they are duplicates or unique, continue to next one
            continue
        else %otherwise add to matrix in the row of the common nodes
        similarity_analysis_y(commonnodes,i) = similarity_analysis_y(commonnodes,i) + 1;
        end

    end
end


%% Old output patterns

old_nodes_fire_normal = readmatrix('Old_OP_normal.csv');

nodes_fired_o = readmatrix('Old_nodes_fired.csv'); %get the output patterns for young learning
nodes_fired_sorted_o = sort(nodes_fired_o); %put them in order

similar6_o = zeros(1,length(nodes_fired_sorted_o(1,:))); %create a matrix to store how many patterns were duplicates

for ii = 1:length(nodes_fired_sorted_o(1,:)) %for all the output patterns

    %verify if the each node firing is the same 
condition1 = find(nodes_fired_sorted_o(1,:) == nodes_fired_sorted_o(1,ii));
condition2 = find(nodes_fired_sorted_o(2,:) == nodes_fired_sorted_o(2,ii));
condition3 = find(nodes_fired_sorted_o(3,:) == nodes_fired_sorted_o(3,ii));
condition4 = find(nodes_fired_sorted_o(4,:) == nodes_fired_sorted_o(4,ii));
condition5 = find(nodes_fired_sorted_o(5,:) == nodes_fired_sorted_o(5,ii));
condition6 = find(nodes_fired_sorted_o(6,:) == nodes_fired_sorted_o(6,ii));
    
    %find whether they overlap across nodes
intersect1 = intersect(condition1,condition2);
intersect2 = intersect(intersect1,condition3);
intersect3 = intersect(intersect2,condition4);
intersect4 = intersect(intersect3,condition5);
intersect5 = intersect(intersect4,condition6);


    %if all the conditions are met, add the number to the matrix, -1 to
    %account for the patterns comparing itself
similar6_o(ii) = length(intersect5) - 1;

end

length(find(similar6_o))
max(similar6_o)

r = 0;
duplicates_o = zeros;

for o = 1:length(nodes_fired_sorted_o)
    for y = 1:length(nodes_fired_sorted_o)
        
         if o == y
            continue
         end
         
        equality = intersect(nodes_fired_sorted_o(:,o), nodes_fired_sorted_o(:,y)); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 6 %if they are duplicates or unique, continue to next one
                disp(o) %show the index of the old OP
                disp(y) %show the index of the young OP
                r = r + 1;
                duplicates_o(r,1) = o;
                duplicates_o(r,2) = y;
        else %otherwise add to matrix in the row of the common nodes
        continue
        end
    end
    if length(duplicates_o(:,1)) == round(length(find(similar6_o))*0.01)
        break
    end
end

%once I've found the duplicates, I will compare their input patterns

for d = 1:length(duplicates_o(:,1))
    first = valid_input_patterns(:,duplicates_o(d,1));
    second = valid_input_patterns(:,duplicates_o(d,2));
    inters = intersect(find(first), find(second));
    sim_dup_o(d) = length(inters);
end

%Uncomment the next line to save matrix as a file (and not have to run this
%again)
% writematrix(sim_dup_o, 'Old_identical_sim.csv')

%plot the figure
    %in case you've already run the script before and only want the
    %figures, run these two lines
% sim_dup_o = readmatrix('Old_identical_sim.csv');
% sim_dup_y = readmatrix('Young_identical_sim.csv');

figure
subplot(1,2,1)
g1 = histcounts(sim_dup_y); %use the young data 
b = bar(g1,'FaceColor',["#EDB120"]); %plot the graph
s = compose('%.1f%%', g1 / sum(g1) * 100); %show the percentages
s1 = compose('%d',g1);
str = [s s1];
yOffset = 1; %position above the bar
yOffset2 = 5;
ycoor = [b.YEndPoints + yOffset, b.YEndPoints + yOffset2];
xcoor = [b.XData, b.XData];
text(xcoor, ycoor, str,'HorizontalAlignment','center','VerticalAlignment','bottom') %alingment of the datatip
set(gca,'XtickLabels',arrayfun(@num2str,[0:5],'UniformOutput',false)); %setting the x-axis
xlabel('Nodes in common',FontSize=15) %print x-axis label
ylabel("Number of patterns",FontSize=15) %print x-axis label
sy1 = sprintf('Sample size = %d', length(sim_dup_y));
text(3, 120, sy1, FontSize=12); %position of writing on the graph
title("Young", FontSize=15)
set(gcf,'color','w');
sgtitle('Common nodes from input patterns with identical output patterns',FontSize=18)

subplot(1,2,2)
g2 = histcounts(sim_dup_o); %use the young data 
b = bar(g2,'FaceColor',["#4DBEEE"]); %plot the graph
so = compose('%.1f%%', g2 / sum(g2) * 100); %show the percentages
so1 = compose('%d',g2);
stro = [so so1];
yOffset = 5; %position above the bar
yOffset2 = 10;
ycoor = [b.YEndPoints + yOffset, b.YEndPoints + yOffset2];
xcoor = [b.XData b.XData];
text(xcoor, ycoor, stro,'HorizontalAlignment','center','VerticalAlignment','bottom') %alingment of the datatip
set(gca,'XtickLabels',arrayfun(@num2str,[0:5],'UniformOutput',false)); %setting the x-axis
xlabel('Nodes in common',FontSize=15) %print x-axis label
ylabel("Number of patterns",FontSize=15) %print x-axis label
s1 = sprintf('Sample size = %d', length(sim_dup_o));
text(3, 120, s1, FontSize=12); %position of writing on the graph
set(gcf,'color','w');
title("Old", FontSize=15)

%% unique patterns

unique_o = zeros(1, length(nodes_fired_sorted_o(1,:)));

for ee = 1:length(nodes_fired_sorted_o(1,:))

    condition1 = find(nodes_fired_sorted_o(1,:) == nodes_fired_sorted_o(1,ee));
    condition2 = find(nodes_fired_sorted_o(2,:) == nodes_fired_sorted_o(2,ee));
    condition3 = find(nodes_fired_sorted_o(3,:) == nodes_fired_sorted_o(3,ee));
    condition4 = find(nodes_fired_sorted_o(4,:) == nodes_fired_sorted_o(4,ee));
    condition5 = find(nodes_fired_sorted_o(5,:) == nodes_fired_sorted_o(5,ee));
    condition6 = find(nodes_fired_sorted_o(6,:) == nodes_fired_sorted_o(6,ee));

    if length(condition1) == 0 && length(condition2) == 0 && length(condition3) == 0 && length(condition4) == 0 && length(condition5) == 0 && length(condition6) == 0
        unique_o(ee) = 1;
    end
end

find(unique_o)

%% Using a representative sample, checking node similarity for 5-1 nodes in common

%now, to determine how many patterns have 5-1 node in common, we must pick
%a representative sample to compare the nodes to.

interest_inputs_o = nodes_fired_sorted_o(:,randomselection); %get the patterns for the randomly selected values 

similarity_analysis_o = zeros(length(similarity_values),length(interest_inputs_o)); 
                                    %matrix with values of interest and
                                    %patterns of interest

for i = 1:length(interest_inputs_o) %random selection
    for ii = 1:length(nodes_fired_sorted_o) %all the patterns
        
        equality = intersect(interest_inputs_o(:,i), nodes_fired_sorted_o(:,ii)); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 0 || commonnodes == 6 %if they are duplicates or unique, continue to next one
            continue
        else %otherwise add to matrix in the row of the common nodes
        similarity_analysis_o(commonnodes,i) = similarity_analysis_o(commonnodes,i) + 1;
        end

    end
end

%% Getting the statistics

max6_o = max(similar6_o);
max6_y = max(similar6new);

%for old learning
for r = 1:length(similarity_analysis_o(:,1))
            my_name = strcat('o',num2str(r)); %give the variable a name
            maximum.(my_name) = max(similarity_analysis_o(r,:)); %assign it to its IP
            minimum.(my_name) = min(similarity_analysis_o(r,:));
            avg.(my_name) = mean(similarity_analysis_o(r,:));
end

%for young learning
for r = 1:length(similarity_analysis_y(:,1))
            my_name = strcat('y',num2str(r)); %give the variable a name
            maximum.(my_name) = max(similarity_analysis_y(r,:)); %assign it to its IP
            minimum.(my_name) = min(similarity_analysis_y(r,:));
            avg.(my_name) = mean(similarity_analysis_y(r,:));
end


