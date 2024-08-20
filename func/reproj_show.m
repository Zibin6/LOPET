function [repro_Image]=reproj_show(K,pose_cur,event_cur,P_3d)

p_pro=K*[pose_cur]*P_3d;
p_pro= [p_pro(1,:)./p_pro(3,:); p_pro(2,:)./p_pro(3,:);ones(1,size(p_pro,2))];
pp1=p_pro(:,1:size(p_pro,2)/2);  pp2=p_pro(:,size(p_pro,2)/2+1:end);
image_cur=zeros(720,1280);
for i=1:size(event_cur,1)
    if event_cur(i,3)==0
        event_cur(i,3)=1;
    end
    if event_cur(i,2)==0
        event_cur(i,2)=1;
    end
    image_cur(event_cur(i,3), event_cur(i,2)) = 255;
end
imshow(image_cur);
hold on
line([pp1(1,:);pp2(1,:)],[pp1(2,:);pp2(2,:)],'color','r','LineWidth',1);
hold off

repro_Image = getframe;
repro_Image = repro_Image.cdata;


end

