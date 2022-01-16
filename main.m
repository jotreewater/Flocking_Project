% CS 455 - Flocking Project

clc,clear
close all

% User Inputs
case_number = 4;
delta_t = .009;
d = 15;
t = 0:delta_t:5;
number_nodes = 100;
% Variables
k_scale = 1.2;
r = k_scale * d;
number_dimensions = 2;
if case_number == 1 || case_number == 2
    area = 50;
    t = 0:delta_t:5;
end
if case_number == 3 || case_number == 4
    area = 150;
    t = 0:delta_t:15;
end
q = area .* rand(number_nodes, number_dimensions);
p = zeros(number_nodes,number_dimensions);
q_old = q;
q_mt = zeros(size(t,2),number_dimensions);
p_mt = zeros(size(t,2),number_dimensions);
q_mt_old = zeros(size(t,2),number_dimensions);
q_mean = zeros(size(t,2),number_dimensions);
p_mean = zeros(size(t,2),number_dimensions);
connectivity = zeros(size(t,2),number_dimensions);
q_all = cell(size(t,2),1);
p_all = cell(size(t,2),1);
n_frames = length(t); % Sets number of frames for the movie
frames = struct('cdata', cell(1, n_frames),'colormap', cell(1,n_frames)); % Preallocates the movie structure
% Iteration Loop
for k = 1:length(t)
    
    if(case_number == 2)
        q_mt(k,:) = [150,150];
        if(k>1)
            q_mt_old(k,:) = q_mt(k-1,:);
            p_mt(k,:) = ((q_mt(k,:) - q_mt_old(k,:)));
        end
        plot(q_mt(:,1),q_mt(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2);
        hold on;
    elseif (case_number == 3)
        q_mt(k,:) = [75 + 50 * t(k),75 + 50 * sin(t(k))];
        if(k>1)
            q_mt_old(k,:) = q_mt(k-1,:);
            p_mt(k,:) = ((q_mt(k,:) - q_mt_old(k,:)));
        end
        plot(q_mt(:,1),q_mt(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2);
        hold on;
    elseif (case_number == 4)
        q_mt(k,:) = [75 + 130*cos(t(k)),75 + 130*sin(t(k))];
        if(k>1)
            q_mt_old(k,:) = q_mt(k-1,:);
            p_mt(k,:) = ((q_mt(k,:) - q_mt_old(k,:)));
        end
        plot(q_mt(:,1),q_mt(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2);
        hold on;
    else
        
    end
    % q,p,u code
    if(case_number == 1)
        [Ui] = Algorithm1(q,p,number_nodes,r,d);
    end
    if(case_number == 2 || case_number == 3 || case_number == 4)
        [Ui] = Algorithm2(q,p,number_nodes,r,d,q_mt,p_mt,k);
    end
    p = ((q - q_old)/delta_t);
    p_all{k} = p;
    q_old = q;
    % Update Positions
    q = q_old + p * delta_t + Ui * delta_t * delta_t / 2;
    [q_ik, A] = Find_Neighbors(q, r, number_nodes);
    if(case_number > 2)
        q_mean(k,:) = mean(q);
        plot(q_mean(:,1),q_mean(:,2),'ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k','MarkerSize',4.2);
        hold on;
    end
    q_all{k} = q;
    connectivity(k) = ((1/number_nodes) * rank(A));
    plot(q(:,1),q(:,2), '.');
    hold on
    for i = 1:number_nodes
        for j = 1:size(q_ik{i})
            plot([q(i,1),q_ik{i}(j,1)],[q(i,2),q_ik{i}(j,2)], 'b');
            hold on;
        end
    end
    plot(q(:,1),q(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
    hold off
    frames(k) = getframe;
end
VW = VideoWriter(sprintf('Case%d.avi', case_number));
open(VW);
writeVideo(VW, frames);
close(VW);
clear frames;

% Post-Simulation Plots

% Trajectory
for i = 1:length(q_all)
    figure(2), plot(q_all{i}(:,1), q_all{i}(:,2), 'k.')
    hold on
end

% Velocity
velocity_graph = zeros(size(t,2),number_nodes);
for i = 1:number_nodes
    for j = 2:length(t)
        temp = norm(p_all{j-1}(i,:) - p_all{j}(i,:));
        velocity_graph(j,i) = temp;
    end
    figure(3),plot(velocity_graph);
    hold on;
end

% Connectivity
figure(4), plot(connectivity);
hold on;

% Functions
function [q_ik, A] = Find_Neighbors(q, r, number_nodes)
q_ik = cell(number_nodes,1);
A = zeros(number_nodes,number_nodes);
for i = 1:number_nodes
    neighbor_number = 1;
    for j = 1:number_nodes
        if((i ~= j) && (norm(q(i,:)-q(j,:)) < r))
            q_ik{i}(neighbor_number,:) = q(j,:);
            neighbor_number = neighbor_number + 1;
            A(i,j) = 1;
        end
    end
end
end
function [Ui] = Algorithm1(q,p,number_nodes,r,d)
Ui = zeros(number_nodes, 2);
c_a1 = 30;
c_a2 = 2 * sqrt(c_a1);
for i = 1:number_nodes
    gradient = 0;
    consensus = 0;
    for j = 1:number_nodes
        gradient = gradient + (Phi_a(Sigma_norm(q(j,:) - q(i,:)), r, d) * N_ij(q(i,:),q(j,:)));
        if(i ~= j)
            a_ij = A_ij(q(i,:),q(j,:),r);
        else
            a_ij = 0;
        end
        consensus = consensus + (a_ij * (p(j,:) - p(i,:)));
    end
    Ui(i,:) = c_a1 * gradient + c_a2 * consensus;
end
end
function phi_a = Phi_a(z,r,d)
r_a = Sigma_norm(r);
d_a = Sigma_norm(d);
phi_a = Bump(z / r_a) * Phi(z - d_a);
end
function phi = Phi(z)
a = 5;
b = 5;
c = abs((a - b)) / sqrt(4 * a * b);
phi = (1/2) * ((a + b) * Sigma_1(z + c) + (a - b));
end
function sigma_1 = Sigma_1(z)
sigma_1 = z / sqrt(1+z^2);
end
function n_ij = N_ij(q_i,q_j)
n_ij = Sigma_gradient(q_j - q_i);
end
function a_ij = A_ij(q_i,q_j,r)
r_a = Sigma_norm(r);
a_ij = Bump(Sigma_norm(q_j - q_i) / r_a);
end
function [Ui] = Algorithm2(q,p,number_nodes,r,d,q_mt,p_mt,k)
Ui = zeros(number_nodes, 2);

c_a1 = 250;
c_a2 = 2 * sqrt(c_a1);
c_mt1 = 20;
c_mt2 = 2 * sqrt(c_mt1);
for i = 1:number_nodes
    gradient = 0;
    consensus = 0;
    for j = 1:number_nodes
        gradient = gradient + (Phi_a(Sigma_norm(q(j,:) - q(i,:)), r, d) * N_ij(q(i,:),q(j,:)));
        if(i ~= j)
            a_ij = A_ij(q(i,:),q(j,:),r);
        else
            a_ij = 0;
        end
        consensus = consensus + (a_ij * (p(j,:) - p(i,:)));
    end
    q_gamma = q(i,:) - q_mt(k,:);
    p_gamma = p(i,:) - p_mt(k,:);
    Ui(i,:) = c_a1 * gradient + c_a2 * consensus - c_mt1 * q_gamma - c_mt2 * p_gamma;
end
end






























