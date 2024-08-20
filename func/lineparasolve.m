function [line] = lineparasolve(point)

  for i=1:size(point,2)/2
 
      line(1:3,i)= cross(point(1:3,i),point(1:3,(size(point,2)/2+i)),1);
      line(1:3,i)= xnorm(line(1:3,i));
  end
end



function V1 = xnorm(V0)

	V1 = V0 .* kron(ones(3,1), 1./sqrt(sum(V0.^2)));
end
