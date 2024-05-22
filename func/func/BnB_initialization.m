function [R_opt,t_opt] = BnB_initialization(point_3d,point_2d,epsilon)
point_2D_a=point_2d(:,1:2:end);
point_2D_b=point_2d(:,2:2:end);

point_3D_a=point_3d(1:3,1:2:end);
point_3D_b=point_3d(1:3,2:2:end);

data_2d_v_non=point_2D_a-point_2D_b;
data_3d_v_non=point_3D_a-point_3D_b;

data_2d_v=data_2d_v_non./vecnorm(data_2d_v_non);
data_3d_v=data_3d_v_non./vecnorm(data_3d_v_non);

data_2d_c=0.5*(point_2D_a+point_2D_b);
data_3d_c=0.5*(point_3D_a+point_3D_b);


data_2d_N_non=cross(data_2d_v,data_2d_c);
data_2d_N=data_2d_N_non./vecnorm(data_2d_N_non);


best_branch=[-pi;-pi;-pi;pi;pi;pi];
M=[best_branch(1:3),0.5*(best_branch(1:3)+best_branch(4:end)),best_branch(4:end)];
for i=1:8
    new_branch(1,i)= M(1,bitget(i,1)+1);
    new_branch(2,i)= M(2,bitget(i,2)+1);
    new_branch(3,i)= M(3,bitget(i,3)+1);
    new_branch(4,i)= M(1,bitget(i,1)+2);
    new_branch(5,i)= M(2,bitget(i,2)+2);
    new_branch(6,i)= M(3,bitget(i,3)+2);
end


parfor j=1:8
[R_opt0{j} e_matrix{j} r_index0{j}]=BnB_for_Rotation(double(data_3d_v),double(data_2d_N),double(epsilon),new_branch(:,j)) ;

end


[rrr r_index]=min(cell2mat(r_index0));
[mindis inlier_label]=min(e_matrix{r_index},[],2);


R_opt=R_opt0{r_index};

data_2d_N2=[];data_3d_c2=[];
for i=1:size(mindis,1)
    if (mindis(i)<=epsilon)
        data_2d_N2=[data_2d_N2 data_2d_N(:,inlier_label(i))];
        data_3d_c2=[data_3d_c2 data_3d_c(:,i)];
    end
end


 [t_opt] = get_trans_outlier(data_2d_N2,data_3d_c2,R_opt);


end

function [modelRANSAC,inlierIdx] = get_robust_trans(A,y)

sampleSize = 3; 
maxDistance = 0.01; 

fitLineFcn = @(data) data(:,1:3)\data(:,4); 
evalLineFcn =@(model, data) (data(:,1:3)*model-data(:,4)).^2;
warning off 'MATLAB:singularMatrix' 
[modelRANSAC, inlierIdx] = ransac([A,y],fitLineFcn,evalLineFcn,sampleSize,maxDistance);
warning on 'MATLAB:singularMatrix'
end




function [t_opt] = get_trans_outlier(data_2d_N,data_3d_c,R_opt)

data_3d_c_rot=R_opt*data_3d_c;
K=-dot(data_2d_N,data_3d_c_rot);
t_opt = get_robust_trans(data_2d_N',K');


end
