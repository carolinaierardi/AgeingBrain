% Script Name: input_pattern)analysis.m
% Date: 16/10/2022
% Author: CMI
% Version: 3.0
% Purpose: we will first conduct an input pattern similarity analysis and
% secondly, we will investigate whether similar input patterns generate
% similar output patterns.
% Notes: there are two parts to this script (second part done for young and
% old learning).

close all
clear all
clc

%% Input pattern analysis

%get matrices needed here from other scripts: 
valid_input_patterns = readmatrix('IP_after_threshold_normal.csv');
valid_output_patterns = readmatrix('OP_after_threshold_normal.csv');
young_nodes_fire_normal = readmatrix('Young_OP_normal.csv');
old_nodes_fire_normal = readmatrix('Old_OP_normal.csv');
over_th_old = readmatrix('Old_nodes_fired.csv');
over_th_young = readmatrix('Young_nodes_fired.csv');

%now, to determine how many patterns have 5-1 node in common, we must pick
%a representative sample to compare the nodes to.

sample = round(length(valid_input_patterns(1,:))*0.01);
%we will select a sample of 1%
randomselection = randi(length(valid_input_patterns(1,:)), [1,sample]);
interest_input = valid_input_patterns(:,randomselection); %get the patterns for the randomly selected values 

similarity_values = [1:5]; %values I want to test for similarity
similarity_analysis_IP = zeros(length(similarity_values),length(interest_input)); 
                                    %matrix with values of interest and
                                    %patterns of interest

for i = 1:length(interest_input) %random selection
    for ii = 1:length(valid_input_patterns) %all the patterns
        
        equality = intersect(find(interest_input(:,i)), find(valid_input_patterns(:,ii))); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 0 || commonnodes == 6 %if they are duplicates or unique, continue to next one
            continue
        else %otherwise add to matrix in the row of the common nodes
        similarity_analysis_IP(commonnodes,i) = similarity_analysis_IP(commonnodes,i) + 1;
        end

    end
end

%getting min, max and mean
for r = 1:length(similarity_analysis_IP(:,1))
            my_name = strcat('n',num2str(r)); %give the variable a name
            maximum.(my_name) = max(similarity_analysis_IP(r,:)); %assign it to its IP
            minimum.(my_name) = min(similarity_analysis_IP(r,:));
            avg.(my_name) = mean(similarity_analysis_IP(r,:));
end


%% Similar input patterns - similar output pattern?
% 
% Step 1) we will search for patterns that appear both in young and old learning
% Step 2) Then, we will look for their IPs and find which other IPs have 5 nodes in
% common. 
% Step 3) Finally, we will get the OPs they generated and test this subset
% for similarity

%Step 1:
% find output patterns common to young and old learning
r = 0;
young_old_OPs = zeros(10,2);

for o = 1:length(over_th_old)
    for y = 1:length(over_th_young)
        
        equality = intersect(over_th_old(:,o), over_th_young(:,y)); %intersection of common nodes
        commonnodes = length(equality); %length of interesection

        if commonnodes == 6 %if they are duplicates or unique, continue to next one
                disp(o) %show the index of the old OP
                disp(y) %show the index of the young OP
                r = r + 1;
                young_old_OPs(r,1) = o;
                young_old_OPs(r,2) = y;
            break
        else %otherwise add to matrix in the row of the common nodes
        continue
        end
    end
    if length(find(young_old_OPs(:,1))) == 10 %we only want the first 10 
        break
    end
end

%Step 2a) now, we get the input patterns that generated them 

o = 0; %initialising indeces
y = 0;
for p = 1:length(young_old_OPs) %for all the found patterns in common
    for yo = 1:length(young_old_OPs(1,:)) %young and old
        if yo == 1 %the first column has old patterns
            o = o + 1; %update index
            my_field = strcat('o',num2str(o)); %give the variable a name
            O_patt.(my_field) = valid_input_patterns(:,young_old_OPs(p,yo)); %assign it to its IP

        else %if it's in the second column
            y = y + 1; %update index
            my_field = strcat('y',num2str(y)); %give variable a name
            Y_patt.(my_field) = valid_input_patterns(:,young_old_OPs(p,yo)); %assign it to its IP

        end
    end
end

%we will only use the first three patterns found
OPIPanalysis_old = [O_patt.o1, O_patt.o2, O_patt.o3];
OPIPanalysis_young = [Y_patt.y1, Y_patt.y2, Y_patt.y3];


%Step 2b) now we want to get the input patterns that have 5 nodes in common for
%these patterns


for a = 1:length(OPIPanalysis_old(1,:))

    x = 1;
    similarity_IPOP_old(:,1,a) = over_th_old(:,young_old_OPs(a,1)); %make the first column the variable we are interested in

    for f = 1:length(valid_input_patterns(1,:)) %I get all the valid input patterns
        
        equality = intersect(find(OPIPanalysis_old(:,a)), find(valid_input_patterns(:,f))); %see how many active nodes are similar
        [similar y] = size(equality); %how many nodes are in common
        
        if similar == 5 %if there are 5 in common, I add to the matrix
            x = x + 1;
            similarity_IPOP_old(:,x,a) = over_th_old(:,f); %add their output pattern
        end
    end
end


%Step 3) now we determine how many nodes in the output patterns are activated

for a = 1:length(OPIPanalysis_old(1,:))
    for ee = 1:length(similarity_IPOP_old(1,:,1))
        for ii = 1:length(similarity_IPOP_old(1,:,1)) %for each of the output patterns

            if ee == ii %if the IPs are the same, don't compare
                n_nodes_old(ee,ii,a) = NaN;
                continue
            elseif similarity_IPOP_old(1,ee,a) == 0 %if there were fewer IPs, don't compare
                n_nodes_old(ee,ii,a) = NaN;
                continue
            elseif similarity_IPOP_old(1,ii,a) == 0 %if there were fewer IPs, don't compare
                n_nodes_old(ee,ii,a) = NaN;
                continue
            end

            equality = intersect(similarity_IPOP_old(:,ee,a), similarity_IPOP_old(:,ii,a)); %compare the output patterns with each output pattern
            n_nodes_old(ee,ii,a) = length(equality);

        end
    end
end


%% For young learning

%Step 2b) now we want to get the input patterns that have 5 nodes in common for
%these patterns

for a = 1:length(OPIPanalysis_young(1,:))

    x = 1;
    similarity_IPOP_young(:,1,a) = over_th_young(:,young_old_OPs(a,2)); %make the first column the variable we are interested in

    for f = 1:length(valid_input_patterns(1,:)) %I get all the valid input patterns
        
        equality = intersect(find(OPIPanalysis_young(:,a)), find(valid_input_patterns(:,f))); %see how many active nodes are similar
        [similar y] = size(equality); %how many nodes are in common
        
        if similar == 5 %if there are 5 in common, I add to the matrix
            x = x + 1;
            similarity_IPOP_young(:,x,a) = over_th_young(:,f); %add their output pattern
        end
    end
end


%Step 3) now we determine how many nodes in the output patterns are activated

for a = 1:length(OPIPanalysis_old(1,:))
    for ee = 1:length(similarity_IPOP_young(1,:,1))
        for ii = 1:length(similarity_IPOP_young(1,:,1)) %for each of the output patterns

            if ee == ii
                n_nodes_young(ee,ii,a) = NaN;
                continue
            elseif similarity_IPOP_young(1,ee,a) == 0
                n_nodes_young(ee,ii,a) = NaN;
                continue
            elseif similarity_IPOP_young(1,ii,a) == 0
                n_nodes_young(ee,ii,a) = NaN;
                continue
            end

            equality = intersect(similarity_IPOP_young(:,ee,a), similarity_IPOP_young(:,ii,a)); %compare the output patterns with each output pattern
            n_nodes_young(ee,ii,a) = length(equality);

        end
    end
end

%% plotting
%initalising loop variables
x = 0;
y = 0;

figure
for graphs = 1:6
    subplot(3,2,graphs)
       % young patterns on the left side and old on the right
    if graphs  == 1 || graphs == 3 || graphs == 5 %young pattern side
        x = x + 1; %keep track of young index
        g1 = histcounts(n_nodes_young(:,:,x)); %use the young data
        b = bar(g1,'FaceColor',["#EDB120"]); %plot the graph
        v = n_nodes_young(1,:,x); %use data of the pattern in question
        s1 = sprintf('IPs - 5 common nodes = %d', length(v(~isnan(v))) ); %how many patterns there were, excluding NA
    else
        y = y + 1; %keep track of old index
        g1 = histcounts(n_nodes_old(:,:,y)); %otherwise, use old data
        b = bar(g1,'FaceColor',["#4DBEEE"]);
        v = n_nodes_old(1,:,y); %use data of the pattern in question
        s1 = sprintf('IPs - 5 common nodes = %d', length(v(~isnan(v)))); %how many patterns there were, excluding NA
    end

    s = compose('%.1f%%', g1 / sum(g1) * 100); %show the percentages
    yOffset = 5; %position above the bar
    text(b.XData, b.YEndPoints + yOffset, s,'HorizontalAlignment','center','VerticalAlignment','bottom') %alingment of the datatip
    set(gca,'XtickLabels',arrayfun(@num2str,[0:5],'UniformOutput',false)); %setting the x-axis
    text(length(similarity_values)/1.3, max(g1)/2, s1); %position of writing on the graph
    if graphs  == 1 || graphs == 3 || graphs == 5
    str = 'Young%d - nodes in common\n'; %x-axis label
        xlabel(sprintf(str,x)) %print x-axis label
    else
    str = 'Old%d - nodes in common\n'; %x-axis label
        xlabel(sprintf(str,y)) %print x-axis label
    end
    ylabel("Output patterns") %y-label is always the same 
    sgtitle("Output pattern similarity for similar input patterns",FontSize=15) %title for graph
    set(gcf,'color','w');
end



