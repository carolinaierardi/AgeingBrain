%Script Name: bipartite_graph_KURF.m
%Date: 29/07/2022
%Author: CMI
%Version: 2.0
%Purpose: create bipartite graph with 30 input and ouput nodes
%Notes: In this version we will investigate the variance and the
%distributions. 

close all
clear all
clc

    %we need a  graph with 30 input and 30 output nodes. Every input node
    %has an edge connecting to every output node. Therefore, each input
    %node has 30 edges flowing from it. The sum of weights of each output
    %node is 1. We need to:
    %a) create ways to assign weights to the input nodes and;
    %b) analyse the variance in the sum of weights reaching the output
    %nodes. 

%% Using uniform distribution

%Creating a bipartite graph 
input_u = zeros(30, 30); %We create a matrix with 2 dimensions representing number of nodes, 
                 % number of edges
                 %values are initiatialised to 0

%in the matrix created above, I need to distribute weights to the edges so 
% that the sum of a row = 1.

[n_node_u, n_edge_u] = size(input_u); %naming the variables that represent the dimensions of the graph


for i = 1:n_node_u %loop for every node in the input nodes
        w_u(i,:) = rand(1,n_edge_u,'double'); %create 30 random numbers for each edge of that node
        w_norm_u = w_u(i,:) / sum(w_u(i,:)); %normalise to sum to 1
        %w_round = round(w_norm_u(i,:), 3, "decimals"); %rounds to 3 decimals
        input_u(i,:) = w_norm_u; %put values into matrix created before
end

%therefore, the sum of each column will be equal to the signal arriving at
%each output node

for e = 1:n_edge_u
output_u(e) = sum(input_u(:,e));
end

disp(output_u); %display the sum of the signal arriving to each output neuron

%calculating some data for this procedure
variance_u = var(output_u);
minimum_u = min(output_u);
maximum_u = max(output_u);

%producing a heatmap of the values found
figure; %set up figure environment
subplot(3,2,3) %subplot grid 1
h_u = heatmap(output_u', 'Title', 'Heatmap - output signal in nodes (Uniform distribution)', FontSize=12);
h_u.ColorScaling = 'scaled';
colormap("winter")
ylabel("Nodes")
%xlabel("Heatmap")
set(gcf,'color','w');

subplot(3,2,1)
iprandu = histogram(input_u, 10, DisplayStyle="bar", FaceColor=[0.2 0.2 0.9]);
text(0.04,250,'Weight assignment distribution curve for edges',"FontSize",16, FontWeight="bold");
%title("Weight assignment distribution curve for edges", "FontSize",18)
ylabel("Number of edges", "FontSize",14)
xlabel("Weight signal", "FontSize",14)
ylim([0 200])
title("Uniform distribution", "FontSize",14)

subplot(3,2,5) %subplot grid 2
bars = histogram(output_u',10, DisplayStyle="bar", FaceColor= [0.2 0.9 0.2]);
%title("Resulting Output signal distribution", "FontSize",18)
ylim([0 10])
text(1.1, 11, "Resulting Output signal distribution", "FontSize",16,"FontWeight","bold")
ylabel("Number of nodes", "FontSize",14)
xlabel("Signal", "FontSize",14)
set(gcf,'color','w');

%% Using normal distribution

%Creating a bipartite graph 
input_n = zeros(30, 30); %We create a matrix with 2 dimensions representing number of nodes, 
                 % number of edges
                 %values are initiatialised to 0

%in the matrix created above, I need to distribute weights to the edges so 
% that the sum of a row = 1.

[n_node_n, n_edge_n] = size(input_n); %naming the variables that represent the dimensions of the graph


for a = 1:n_node_n %loop for every node in the input nodes
        w_n(a,:) = 0.3*randn(1,n_edge_n, 'double') + 1 ; %create 30 random numbers for each edge of that node
        w_norm_n = w_n(a,:) / sum(w_n(a,:)); %normalise to sum to 1
        input_n(a,:) = w_norm_n; %put values into matrix created before
end

%therefore, the sum of each column will be equal to the signal arriving at
%each output node

for b = 1:n_edge_n
output_n(b) = sum(input_n(:,b));
end

disp(output_n); %display the sum of the signal arriving to each output neuron

%calculating some data for this procedure
variance_n = var(output_n);
minimum_n = min(output_n);
maximum_n = max(output_n);

%producing a heatmap of the values found


subplot(3,2,4)
h_n = heatmap(output_n','Title', 'Heatmap - output signal in nodes (Normal distribution)', FontSize=12);
h_n.ColorScaling = 'scaled';
colormap("winter")
ylabel("Nodes")
%xlabel("Heatmap")
set(gcf,'color','w');

subplot(3,2,2)
iprandn = histogram(input_n, 10, DisplayStyle="bar", FaceColor=[0.2 0.2 0.9]);
%title("Weight assignment distribution curve for edges", "FontSize",18)
ylabel("Number of edges", "FontSize",14)
xlabel("Weight signal", "FontSize",14)
ylim([0 300])
title("Normal distribution", "FontSize",14)

subplot(3,2,6)
bars = histogram(output_n',10, DisplayStyle="bar", FaceColor= [0.2 0.9 0.2]);
%title("Resulting Output signal distribution", "FontSize",18)
ylim([0 10])
ylabel("Number of nodes", "FontSize",14)
xlabel("Signal", "FontSize",14)
set(gcf,'color','w');


%%Uncomment the next two lines to generate new edge weights to use in next scripts
%writematrix(input_u,"Inputs_uniform_dis.csv")
%writematrix(input_n, "Inputs_normal_dis.csv")

                %%% END OF SCRIPT %%%
