%%  The demo code for Line-based 6-DoF Object Pose Estimation and Tracking With an Event Camera
clc
clear
close all
addpath(genpath(pwd))
load('line_parameter.mat')

%% Pose Initialization

p = gcp('nocreate');
    if isempty(p) && license('test', 'Distrib_Computing_Toolbox') 
        % Parallel Computing Toolbox is licensed, not running locally, and no pool is open  
        disp('Starting parallel pool...');  
        Mypar = parpool; % Start a parallel pool  
    else  
        % Either a pool is already open or parallel computing is not supported  
        if ~isempty(p)  
            disp('Parallel pool is already open.');  
        else  
            disp('Parallel computing is not supported or cannot be started.');  
        end  
    end  


fprintf('Start pose initialization...\n');

tic
[R_opt,t_opt_] =BnB_initialization(X_W,x_c,0.0299);
time_initial=toc
pose_cur=[R_opt,t_opt_];
% delete(Mypar);% Parallel pool closed.

%% Pose Tracking
fprintf('Start pose tracking...\n');
j=1;

for i=1:200

    event_window_name=['.\events\',num2str(i),'.mat'];
    load(event_window_name);

    p_cur=K*pose_cur*X_W;
    p_cur= [p_cur(1,:)./p_cur(3,:); (p_cur(2,:)./p_cur(3,:));ones(1,size(p_cur,2))];
    [event_cluster] = point_cluster(p_cur,event_cur,10);

    tic
    [pose_cur] =pose_optim(event_cluster, pose_cur(1:3,1:3),pose_cur(1:3,4),K,X_W);
    time_m(j,1)=toc;

    [repro_Image]=reproj_show(K,pose_cur,event_cur0,X_W);
    imshow(repro_Image);
    j=j+1;

end

%% Time-consuming
fprintf('The time consumed for pose initialization is %d s\n', time_initial);

time_optimization=median(time_m);
fprintf('The time consumed for pose optimization is %d s\n', time_optimization);
