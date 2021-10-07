% Joseph Trierweiler
% CS455 - Mobile Sensor Networks
% Project 1

clc,clear
close all

% User Defined Inputs

p_case = 1;
d = 15; % Desired Distance between nodes
epsilon = 0.1; % Determines sigma norm strength (between 0 and 1)
n_nodes = 50; % Desired number of sensor nodes
n_dimensions = 2; % Desired number of dimensions
d_time = 0.1; %Set time step
time = 0:d_time:7.5;% Set simulation time

% Static Variables

k = 1.2; % Scaling Factor
r = k * d; % Node Interaction Range

% Variables
if (p_case == 1 || p_case == 2)
    area = 50; % Initial deployment range of nodes (area * area)
else
    area = 150; % Initial deployment range of nodes (area * area)
end
nodes = area.*rand(n_nodes, n_dimensions) + area.*repmat([0 1], n_nodes, 1); % Generates Inital nodes
v_nodes = zeros(n_nodes,n_dimensions); % Node velocity matrix
p_nodes = nodes; % Previous node coordinates
p_a_position = zeros(size(time,2), n_dimensions); % Previous average posotion of nodes
p_a_velocity = zeros(size(time,2), n_dimensions); % Previous average velocity of nodes
connectivity = (1000); % How connected the MSN is
p_positions = cell(size(time,2), n_nodes); % Cell array for recording node positions
p_velocities = cell(size(time,2), n_nodes); % Cell array for recording node velocities
p_t_nodes = cell(size(time,2), n_nodes); % Cell array for recording previous taraget nodes
t_n_velocity = [0 0]; % Velocity of the target node
t_node = [size(time,2), n_nodes];

n_frames = length(time); % Sets number of frames for the movie
frames = struct('cdata', cell(1, n_frames),'colormap', cell(1,n_frames)); % Preallocates the movie structure

% Main loop
for i = 1:length(time)
    % Target Code
    if(p_case > 2)
        if (p_case == 3)
            % Move target node in a sin wave
            t_node(i,:) = [50 + 50*time(i),295 - 50*sin(time(i))];
            % Calculate the velocity of the target node
            if(i > 1)
                p_t_nodes {i} = t_node;
                t_n_velocity(i,:) = (t_node(i,:) - t_node(i-1,:)) / d_time;
            end
        elseif (p_case == 4)
            % Move the target node in a circle
            t_node(i,:) = [310 - 160*cos(time(i)),255 + 160*sin(time(i))];
            % Calculate the velocity of the target node
            if(i > 1)
                p_t_nodes {i} = t_node;
                t_n_velocity(i,:) = (t_node(i,:) - t_node(i-1,:)) / d_time;
            end
        end
        % Plot Target Node
        if(i > 1)
            figure(1),plot(p_t_nodes{i}(:,1),p_t_nodes{i}(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2)
            hold on
        end
        
    elseif(p_case == 2)
        t_node = [150 150]; % Position of the target node
        figure(1),plot(t_node(:,1),t_node(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2)
        hold on
    end
    % Node Code (haha)
    if p_case == 1
        [a_neighbors, a_matrix] = Find_A_Neighbors(nodes, r);
        [n_accelerations] = Algorithm1(n_nodes, nodes, a_neighbors, n_dimensions, epsilon, r, d, v_nodes);
    elseif p_case == 2
        [a_neighbors, a_matrix] = Find_A_Neighbors(nodes, r);
        [n_accelerations] = Algorithm2(p_case, n_nodes, nodes, a_neighbors, n_dimensions, epsilon, r, d, t_node, t_n_velocity, v_nodes);
    elseif p_case == 3 || p_case == 4
        [a_neighbors, a_matrix] = Find_A_Neighbors(nodes, r);
        [n_accelerations] = Algorithm2(p_case, n_nodes, nodes, a_neighbors, n_dimensions, epsilon, r, d, t_node(i,:), t_n_velocity(i,:), v_nodes);
    end
    
    % Record Changes
    v_nodes = (nodes - p_nodes) / d_time; % Calculate Velocity
    p_velocities{i} = v_nodes;  % Record Node Velocity
    p_nodes = nodes; % Record Node Position
    nodes = p_nodes + v_nodes*d_time + .5 * n_accelerations * d_time * d_time;
    p_a_position(i,:) = mean(nodes); % Record Average Node Position
    p_a_velocity(i,:) = mean(v_nodes); % Record Average Node Acceleration
    p_positions{i} = nodes; % Record Node Positions into Cell Array
    connectivity(i)= (1 / n_nodes) * rank(a_matrix); % Record Node Connectivity
    
    if (p_case > 1)
        % Plot Average Node Position
        hold on
        plot(p_a_position(:,1),p_a_position(:,2),'ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k','MarkerSize',4.2)
        
    end
    
    % Plot Nodes
    plot(nodes(:,1),nodes(:,2), '.')
    hold off
    % plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    
    % Draw Lines
    for j = 1:n_nodes
        temp_node = nodes(a_neighbors{j},:);
        for k = 1:size(nodes(a_neighbors{j},1))
            line([nodes(j,1),temp_node(k,1)],[nodes(j,2),temp_node(k,2)])
        end
    end
    % Record Figure Frame for movie
    frames(i) = getframe;
end

% Movie Code

VW = VideoWriter(sprintf('Case%d.avi', p_case));
open(VW);
writeVideo(VW, frames);
close(VW);
clear frames;

% Trajectory Plot
for i = 2:length(p_positions)
    hold on
    figure(2), plot(p_positions{i}(:,1), p_positions{i}(:,2), 'k.')
end

% Velocity Plot
p_each_nodes = (length(time));
for i = 2:size(time,2)
    for j = 1:n_nodes
        if (j == 1)
            p_each_nodes(i) = norm(p_velocities{i}(j,:));
            hold on
            figure(3), plot(p_each_nodes, 'b')
        end
    end
end

% Connectivity plot
hold on
figure(4), plot(connectivity)

if (p_case == 3)
    figure(5), plot(p_a_position(:,1), p_a_position(:,2),'k.')
    hold on
    plot(t_node(:,1), t_node(:,2), 'r.')
end
function signorm = sigma_norm(input)
epsilon = .1;
signorm = (1 / epsilon) * (sqrt(1 + epsilon * (norm(input,2)^2)) - 1);
end
function out = bump(z)
if (z >= 0) && (z < .2)
    out = 1;
elseif (z >= .2) && (z <= 1)
    out = .5*(1 + cos(pi*((z - .2) / .8)));
else
    out = 0;
end
end
function [a_neighbors, a_matrix] = Find_A_Neighbors(nodes, r)
n_nodes = size(nodes,1);
dif = cell(n_nodes,1);
distance_alpha = zeros(n_nodes,n_nodes);
temp_d = zeros(n_nodes,1);
a_neighbors = cell(n_nodes,1);

for i = 1:n_nodes
    dif{i} = repmat(nodes(i,:),n_nodes,1) - nodes;
    temp_node = dif{i};
    
    for j = 1:n_nodes
        temp_d(j,:) = norm(temp_node(j,:));
    end
    
    distance_alpha(i,:)= temp_d;
end

for k = 1:n_nodes
    a_neighbors{k} = find(distance_alpha(:,k) < r & distance_alpha(:,k) ~= 0); %find the neighbors of agent i
end

a_matrix = zeros(n_nodes,n_nodes);

for i = 1:n_nodes
    for j = 1:n_nodes
        if i ~= j
            dist_2nodes = norm(nodes(j,:) - nodes(i,:));
            if dist_2nodes < r && dist_2nodes ~= 0
                a_matrix(i,j) = 1;
            end
        end
    end
end
end
function [n_accelerations] = Algorithm1(n_nodes, nodes, a_neighbors, n_dimensions, epsilon, r, d, v_nodes)
ca1 = 30;
ca2 = 2*sqrt(ca1);
a = 5;
b = 5;
c = abs(a - b) / sqrt(4*a*b);
r_sig = sigma_norm(r);
d_sig = sigma_norm(d);

n_ij = zeros(n_nodes,n_nodes,n_dimensions);
for i = 1:n_nodes
    for j = 1:n_nodes
        q = norm(nodes(j,:) - nodes(i,:));
        sig_grad = (nodes(j,:) - nodes(i,:)) / (1 + epsilon * sigma_norm(nodes(j,:) - nodes(i,:)));
        
        if q < r && q ~= 0 
            n_ij(i,j,:) = sig_grad;
        end
    end
end

U = zeros(n_nodes, n_dimensions);

gradient = 0;
conscensus = 0;
a_ij = zeros(n_nodes,n_nodes);
for i = 1:n_nodes
    for j = 1:size(a_neighbors{i})
        Nei_val = a_neighbors{i}(j);
        if(i ~= Nei_val)
            z = sigma_norm(nodes(Nei_val,:) - nodes(i,:));
            z_phi = z - d_sig; 
            rho_h = bump(z / r_sig);
            sigmoid = (z_phi + c) / sqrt(1 + (z_phi + c)^2);
            phi = .5 * ((a + b) * sigmoid + (a - b));
            
            phi_alpha = rho_h * phi;
            
            a_ij(i,Nei_val) = rho_h;
            
            gradient = phi_alpha * [n_ij(i,Nei_val,1) n_ij(i,Nei_val,2)];
            conscensus = a_ij(i,Nei_val) * (v_nodes(Nei_val,:) - v_nodes(i,:));
        end
    end
    
    U(i,:) = (ca1 * gradient) + (ca2 * conscensus);
end

n_accelerations = U;
end
function [n_accelerations] = Algorithm2(p_case, n_nodes, nodes, a_neighbors, n_dimensions, epsilon, r, d, t_node, t_n_velocity, v_nodes)
ca1 = 30;
ca2 = 2*sqrt(ca1);
c_mt1 = 1.1;
c_mt2 = 2*sqrt(c_mt1);
a = 5;
b = 5;
c = abs(a - b) / sqrt(4*a*b);
r_sig = sigma_norm(r);
d_sig = sigma_norm(d);
% phi_a

% phi

% n_ij component
n_ij = zeros(n_nodes,n_nodes,n_dimensions); 
for i = 1:n_nodes
    for j = 1:n_nodes
        q = norm(nodes(j,:) - nodes(i,:));
        sig_grad = (nodes(j,:) - nodes(i,:)) / ((1 + epsilon * sqrt(sigma_norm(nodes(j,:) - nodes(i,:))^2)));
        
        if (q < r) && (q ~= 0) 
            n_ij(i,j,:) = sig_grad;
        end
    end
end

U = zeros(n_nodes, n_dimensions);

conscensus = 0;
a_ij = zeros(n_nodes,n_nodes);
for i = 1:n_nodes 
    gradient = 0;
    for j = 1:size(a_neighbors{i})
        Nei_val = a_neighbors{i}(j);
        if(i ~= Nei_val)
            z = sigma_norm(nodes(Nei_val,:) - nodes(i,:));
            z_phi = z - d_sig;
            phi_bump = bump(z / r_sig);
            sigmoid = (z_phi + c) / sqrt(1 + (z_phi + c)^2);
            phi = .5 * ((a + b) * sigmoid + (a - b));
            
            phi_alpha = phi_bump * phi;
            
            a_ij(i,Nei_val) = phi_bump;
            
            gradient = phi_alpha * [n_ij(i,Nei_val,1) n_ij(i,Nei_val,2)];
            conscensus = a_ij(i,Nei_val) * (v_nodes(Nei_val,:) - v_nodes(i,:));
        end
    end
    
    p = 0;
    if p_case ~= 2
        p = -1 * (c_mt2) * (v_nodes(i,:) - t_n_velocity);
    end
    
    g_force = -1 * (c_mt1) * (nodes(i,:) - t_node) + p;
    a_force = (ca1 * gradient) + (ca2 * conscensus);
    
    U(i,:) = a_force + g_force;
end


n_accelerations = U;
end





