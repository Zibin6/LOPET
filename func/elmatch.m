function [event_cluster event_cur] = elmatch(K,pose_cur,k,cube)

P_3d=[0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000
0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000
-0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	-0.125000000000000	0.125000000000000	-0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000	0.125000000000000
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1];

    p_cur=K*pose_cur*P_3d;
    p_cur= [p_cur(1,:)./p_cur(3,:); (p_cur(2,:)./p_cur(3,:));ones(1,size(p_cur,2))];

    [~, index] = min(abs(k - cube(:,1)));
    start_index = index(1) - 10000;
    end_index = index(1) + 10000;

    event_cur=cube(start_index:end_index,:);


    [event_cluster] = point_cluster(p_cur,event_cur,10);
end

function [event_cluster] = point_cluster(p_cur,event,pix_threshold)
   event_cluster=[];

   lines_cur=lineparasolve(p_cur);

   part_num=size(lines_cur,2);

   for i=1:size(event,1)


       for j=1:size(lines_cur,2)

         point_distance(j)=abs(event(i,2)*lines_cur(1,j) + event(i,3)*lines_cur(2,j) + lines_cur(3,j))/ sqrt(lines_cur(1,j)^2+lines_cur(2,j)^2);

       end
           
        [min_dis,min_row] = min(point_distance);

        if min_dis < pix_threshold  
            
            point_mid=(p_cur(1:2,min_row)+p_cur(1:2,(part_num+min_row)))./2;
            distance_end=norm(p_cur(1:2,min_row)-p_cur(1:2,(part_num+min_row)));
            point_cur=[event(i,2);event(i,3)];
            distance_p_mid=norm(point_cur-point_mid);

            if distance_p_mid<0.55*distance_end 
              event_new     = [event(i,2);event(i,3);1;min_row];
              event_cluster = [event_cluster,event_new];
            else 
              continue
            end
        else
            continue
        end
   end
end


function [line] = lineparasolve(point)

  for i=1:size(point,2)/2
 
      line(1:3,i)= cross(point(1:3,i),point(1:3,(size(point,2)/2+i)),1);
      line(1:3,i)= xnorm(line(1:3,i));
  end
end



function V1 = xnorm(V0)

	V1 = V0 .* kron(ones(3,1), 1./sqrt(sum(V0.^2)));
end
