%Script Name: ip_generation_threshold3.m
%Date: 12/08/2022
%Author: CMI
%Version: 3.0
%Purpose: generate unique input patterns and determine appropriate
%threshold - learning processes
%Notes: THIS IS THE GOOD SCRIPT!!!
    % weighting update adds up between iterations, old learning uses
    % anatomical proximity for selection of the second node
    % We need the RouletteWheelSelection.m file to run this.

close all
clear all
clc

input_u = readmatrix('Inputs_normal_dis.csv'); %weight assignment using normal distribution

%% 3. Input patterns
total_nodes = 30; %our model determines there are 30 total nodes
active_nodes = 6; %our model determines there are 6 nodes firing

p = nchoosek(1:total_nodes, active_nodes); %all combinations of 6 in 30 to determine input patterns
ip = zeros(total_nodes, length(p(:,1))); %create a matrix of IPs of 0s

for i = 1:length(ip(1,:))

ip(p(i,:),i) = 1; %transform active nodes into 1, meaning they fire

end

%here, we have generated a matrix of 30 rows, representing each node and
%593775 columns, representing the input patterns. When a node ==1, it fires
%and when == 0, it does not fire. 
%Now, we need to associate the values the nodes that fire with the weights
%they put in each output node and the sum of signal reaching the 30 output
%nodes in that case. 

newinput = zeros(total_nodes,total_nodes,length(p(:,1))); %I created a matrix of 3 dimensions, with each explained below

[nodes, edges, patterns] = size(newinput); %assigning names to the variables

    %assign the weights to the input patterns 

for h = 1:patterns %for each input pattern 
    for e = 1:nodes %go into each node,...
        if ip(e,h) == 1 %if it fires for that specific pattern ...
            newinput(e,:,h) = input_u(e,:); %assign the weight attributed before
        end
    end
end

    %determine the signal reaching the output nodes for an input pattern

for d = 1:patterns %for each pattern
    for c = 1:nodes % and each node
        newoutput(c,d) = sum(newinput(:,c,d)); %sum the signal reaching each output node
    end
end

%% 4. Determining the threshold

% %As of now, we have 593775 input patterns. We want to determine a threshold
% %so that only 125.000 - 150.000 of these patterns actually reach the output
% %nodes.

th_values = [0.2, 0.2405, 0.28, 0.275, 0.278];% attempt to find a threshold

th_mat = repmat(max(newoutput),length(th_values),1); %get the maximum signal for each IP
valid_th = zeros(length(th_values),1);
 % attempt to find a threshold
 for t = 1:length(th_values)
     for b = 1:patterns %for each input pattern

         if th_mat(t, b) <= th_values(t) %if the maximum value is below the threshold,
             th_mat(t, b) = true; %it is a valid pattern
         else
             th_mat(t, b) = false; %otherwise not valid
         end
     end
     sz = size(find(th_mat(t,:))); %this will give the valid IPs
     valid_th(t) = sz(2); %we store them here to see which threshold worked best
 end

chosen_th = th_values(2); %the second threshold value worked best
valid_ip_index = find(th_mat(2,:)); %index of valid input patterns
                                                    
%plotting the thresholds
graph = max(newoutput); %I get the maximum signal reaching each of the nodes for a given IP

figure 
subplot(1,2,1)
histogram(graph', FaceColor=[0.3 0.9 0.9]) %plot the trans matrix
xlim([0.2 0.32])
xline(chosen_th, Color='red',LineWidth=3) %to show the threshold determined 
title("Maximum signal reaching output nodes in all input patterns", FontSize=18)
xlabel("Maximum signal", FontSize=15)
ylabel("Input Patterns", FontSize=15, FontName='Arial')
s1 = sprintf('Threshold = %.4f',chosen_th);
text(0.28, 2000, s1, "FontSize",13)
set(gcf,'color','w');

%find only the valid input patterns, those with a maximum signal below the
%threshold
valid_input_patterns = ip(:,valid_ip_index);
valid_output = newoutput(:, valid_ip_index);
valid_newinput = newinput(:,:,valid_ip_index);

% getting statistics for the values in the output patterns
variance_valid = max(valid_output) - min(valid_output);
avg_var = mean(variance_valid);
sd_var = std(variance_valid);

%plot the variance in output patterns
subplot(1,2,2)
histogram(variance_valid, FaceColor = [0.1 0.9 0.3])
title("Output variance across each valid input pattern", FontSize=18)
xlabel("Signal variance in output nodes",FontSize=15)
ylabel("Input Patterns", FontSize=15)
s1 = sprintf('Mean variance = %.4f',avg_var); %how many patterns there were, excluding NA
s2 = sprintf('SD variance = %.4f',sd_var);
text(0.11, 2000, s1, "FontSize",13); %position of writing on the graph
text(0.11, 1900, s2, "FontSize",13)
set(gcf,'color','w');

%store the matrices in .csv files
%Uncomment the next three lines to store new values for the IPs for next
%script
% writematrix(valid_newinput, "IPs_after_threshold_weights_normal.csv")
% writematrix(valid_input_patterns, "IP_after_threshold_normal.csv")
% writematrix(valid_output, "OP_after_threshold_normal.csv")

%% Learning processes

    %Doing the young learning process. During this, the weight of an edge
    %is updated according to a function ew1 +(1-ew2)/2, where ew1 is the
    %incoming edge and ew2 is another randomly selected incoming edge. 

    %We wish to determine how many steps does it take to reach the
    %threshold of 0.282 for old and young learning. For that, we will
    %first randomly select one of the 6 activated edges and subsequently
    %select one of the 30 edges 
    
[value, order] = sort(valid_output, 'descend'); 
%the order output of this gives me a matrix with the output nodes' signal sorted from highest to lowest. 
%if the node with the highest signal is the 13th one, it will appear in the
%first row. 


young_n_iterations = zeros(1,length(order(1,:))); %this is to store how many iterations it took to complete learning
over_th_y = zeros(active_nodes,length(order(1,:))); %to store the nodes that fired 
young_nodes_fired = zeros(active_nodes,length(order(1,:)));

for a = 1:length(order(1,:)) %for each input pattern

    loopcount = 0; %iinitialising the amount of iterations
    v = 0;
    for i = 1:nodes %for all the nodes
    chosen_node = RouletteWheelSelection(valid_output(:,a));
    if length(find(young_nodes_fired(:,a) == chosen_node)) ~= 0 %if node has been previously selected, choose another one
        continue
    else
        v = v + 1;
        young_nodes_fired(v,a) = chosen_node; %otherwise, store node for the record

    end

        all = valid_output(chosen_node,a); %initialising the updated weight

        activated = find(valid_newinput(:,chosen_node,a)); %find the active nodes in an input pattern
        all_active = valid_newinput(activated, chosen_node, a); %put them in ascending order

        chosen_edge = RouletteWheelSelection(all_active);
        chosen_edge = all_active(chosen_edge);

        weight = chosen_edge + (1-chosen_edge)/2;

        all = valid_output(chosen_node,a) + weight; %we add this weight to the signal in the output node already
        % if this is >= than the threshold established, we end here.
        % Otherwise, we continue to the next highest node.

        loopcount = loopcount + 1; %counting iterations.

        if all >= chosen_th
            over_th_y(v,a) = all; % if the value is over the threshold, add to the matrix
        end

        if length(find(over_th_y(:,a))) == 6 %if the last value is not 0, stop the loop, go to next input pattern
        break
        end
    
end
    young_n_iterations(a) = loopcount;
end

%The following loop stores output patterns in a vector of 0s and 1s,
%uncomment only if necessary.
%young_learning_firing = zeros(length(active_nodes), length(young_n_iterations(1,:))); %output patterns
% for n = 1:length(young_learning_firing(1,:))
% young_learning_firing(young_nodes_fired(:,n),n) = 1; %transform them into 1
% end

%store the values in a .csv file
%Uncomment next two lines to store values in a .csv for next script
% writematrix(young_learning_firing, 'Young_OP_normal.csv');
% writematrix(young_nodes_fired, 'Young_nodes_fired.csv')

%% For old learning

 %we will choose the first node through roulette selection and the
 %second by anatomical proximity. 
 %We will determine the direction on the web of nodes via a coin flipping simulation

coin = 0.5; %for anatomical selection of the second edge in a node
old_n_iterations = zeros(1,length(order(1,:))); %to store how many iterations learning took
excluded_nodes = zeros(active_nodes,(length(order(1,:)))); %to store which nodes were excluded
over_th_o = zeros(active_nodes,length(order(1,:))); %to store the nodes that fired 
old_nodes_fired = zeros(active_nodes,length(order(1,:))); %matrix with the index of the nodes fired

for a = 1:length(order) %for each input pattern

    loopcount = 0; % amount of iterations

    %CHOOSE A NODE - ROULETTE SELECTION
    for i = 1:nodes %for all the nodes
        chosen_node = RouletteWheelSelection(valid_output(:,a));
        if length(find(old_nodes_fired(:,a) == chosen_node)) ~= 0 %if node has been previously selected, choose another one
            continue
        else
            old_nodes_fired(i,a) = chosen_node; %otherwise, store node for the record
        end


        activated = find(valid_newinput(:,chosen_node,a));%get the active nodes the IP
        all_active = valid_newinput(activated,chosen_node,a); %get the values for the active node

                        %CHOOSE THE FIRST EDGE - ROULETTE SELECTION
        
        chosen_edge = RouletteWheelSelection(all_active);
        chosen_edge = all_active(chosen_edge);
        all = valid_output(chosen_node,a);
       
                        %CHOOSE SECOND EDGE - ANATOMICAL PROXIMITY
 

       %initialising values
        x = rand; %choose a random number between 1 and 0 to define the direction we will go
        s = 0; %to adjust for the position of the selection if threshold isn't met after an iteration
        second_edges = zeros;
        while all <= chosen_th %we will choose nodes until the threshold is met
            all = valid_output(chosen_node,a) + chosen_edge;
            if x >= coin
                %all = valid_output(node_order(i),a) + chosen_edge;
                all_active_circ = [all_active; all_active]; %treat the array as circular
                chosen_loc = find(chosen_edge == all_active); %where does the chosen node match the position in IP
                s = s + 1; %update the positioning

                if s == 6 %in case all edges are chosen and the threshold isn't met, exclude node
                    loopcount = loopcount - 5;
                    excluded_nodes(i,a) = chosen_node;
                    break
                end

                second_edge = all_active_circ(chosen_loc + s); %if greater than 0.5, I go down on the vector
                second_edges(s) = second_edge; %it put all the edges that have been selected
                update = sum(second_edges); %and add them to produce the iterations over the loops
                all = all + update;
            else
                inverse_all_ac = flipud(all_active); %invert the order of active nodes
                inverse_all_circ = [inverse_all_ac; inverse_all_ac]; %treat the array as circular
                chosen_loc = find(chosen_edge == inverse_all_ac); %find the location of the chosen edge
                s = s + 1;
                if s == 6
                    loopcount = loopcount - 5;
                    excluded_nodes(i,a) = chosen_node;
                    break
                end
                second_edge = inverse_all_circ(chosen_loc + s); %use the inverse vector to go down the order
                second_edges(s) = second_edge;
                update = sum(second_edges);
                all = all + update;
            end

            loopcount = loopcount + 1; %counting iterations.

            if all >= chosen_th
                over_th_o(i,a) = chosen_node; % if the value is over the threshold, add to the matrix
            end
        end

        if length(find(over_th_o(:,a))) == 6 %if the last value is not 0, stop the loop, go to next input pattern
            break
        end
    end

    old_n_iterations(a) = loopcount; %to store the iterations for each pattern in a matrix

end

%to later measure old output pattern similarity analysis
old_learning_firing = zeros(nodes, length(old_n_iterations(1,:))); %creating a matrix to see the nodes that fire
over_th_old = zeros(active_nodes, length(old_n_iterations(1,:))); %putting in the index of the nodes


for m = 1:length(old_nodes_fired(1,:))
    over_th_old(:,m) = old_nodes_fired(find(over_th_o(:,m)),m); %only have the fired nodes in the matrix
    %The next loop creates a vector of 0s and 1s with the output patterns,
    %only uncomment if necessary. 
%     for n = 1:length(old_nodes_fired(:,1))
%         old_learning_firing(old_nodes_fired(n,m),m) = 1; %tranform into a vector of 1s and 0s
%     end
end

%store values in .csv files
%Uncomment next two lines to run the next script
% writematrix(old_learning_firing, 'Old_OP_normal.csv')
% writematrix(over_th_old, 'Old_nodes_fired.csv')

%% Figure generation - learning processes
%plotting the number of iterations for old and young learning
figure
younggraph = histcounts(young_n_iterations); %bin counts for the graph
b = bar(younggraph,'FaceColor',["#EDB120"]); %plotting the graph
s = compose('%.1f%%', younggraph / sum(younggraph) * 100); %show % at the datatips
s1 = compose('%d',younggraph); %show actual number at the datatips
yOffset = 5; %position above the end of the bar where numbering starts(1)
yOffset2 = 5000; %position above the end of the bar where numbering starts(2)
str = [s s1]; %string with both datatip info
xcoor = [b.XData b.XData]; %x coordinates for datatips
ycoor = [b.YEndPoints + yOffset, b.YEndPoints + yOffset2]; %y coordinates for datatips
text(xcoor, ycoor, str,'HorizontalAlignment','center','VerticalAlignment','bottom') %position and line for the tips
set(gca,'XtickLabels',arrayfun(@num2str,[6:15],'UniformOutput',false)); %set the x-axis values

hold on %we want old learning on the same graph
oldg = histcounts(old_n_iterations); %bin counts
b = bar(oldg,'FaceColor',["#4DBEEE"]); %plotting graph
s = compose('%.1f%%', oldg / sum(oldg) * 100); %show percentages
s1 = compose('%d',oldg); %show actual number
str = [s s1]; %string with both info
xcoor = [b.XData b.XData]; %x coordinates for text
ycoor = [b.YEndPoints + yOffset, b.YEndPoints + yOffset2]; %y-coordinates for text
text(xcoor, ycoor, str,'HorizontalAlignment','center','VerticalAlignment','bottom') %position for text
set(gca,'XTick',1:length(oldg),'XtickLabels',arrayfun(@num2str,[6:15],'UniformOutput',false)); %setting x-axis
xlabel("Number of iterations","FontSize",15) %x-label
ylabel("Valid input patterns", "FontSize",15) %y-label
ylim([0 150000]) %limit for y-axis, better visualisation
legend("Young", "Old", "FontSize",13) %legends for the different kinds of learning
title("Number of Iterations for Young and Old learning",FontSize=18) %title for the graph
set(gcf,'color','w');

                              %%% END OF SCRIPT%%%
