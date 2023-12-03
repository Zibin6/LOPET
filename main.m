
clc
clear
close all
addpath("func\")
load("demo.mat");


%Camera intrinsic matrix
K=[1654 0 651;0 1653 368;0 0 1];
    
%Initial pose estimation

fprintf('Start initial pose estimation. Please wait.');

pose_cur=initial_mypose();

fprintf('\nThe initial pose estimation is completed, and the result is:\n');

pose_cur

fprintf('\nStart pose tracking');


 for t=cube(11000,1):0.01:cube(end,1)

     if t<(cube(end,1)-0.01)

         % Event - line matching

         [event_cluster event_cur]=elmatch(K,pose_cur,t,cube);
 
         %Pose optimization

         fprintf('\nThe optimized pose is:');    

         [pose_cur] =pose_optim(event_cluster, pose_cur(1:3,1:3),pose_cur(1:3,4),K)

         %Pose reprojection

         [repro_Image]=reproj_show(K,pose_cur,event_cur);

         imshow(repro_Image); 

     else

         fprintf('\nDue to the size restrictions for file uploads. The demo ends here.');

         return;

     end
  
end
