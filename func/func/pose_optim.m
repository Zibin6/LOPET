function [ R_final] =pose_opti( pt, R,T,K,P);


mx = 6;
x_init = zeros(mx, 1);
ang_init = rodrigues(R);
x_init(1:3) = ang_init;
x_init(4:6) = T(1:3,1);
options = optimset;
options = optimoptions("lsqnonlin","Display","none");
options.Algorithm = 'levenberg-marquardt';
options.MaxFunEvals = 1000;
options.TolFun = 1e-3;
options.TolX = 1e-3;
options.MaxIter = 50;
[x_optim,resnorm,residual] =lsqnonlin(@(x) ObjFunReprojErrLine( x, pt, P,K), x_init, [], [], options);
R_optim = rodrigues([x_optim(1) x_optim(2) x_optim(3)]);
T_optim = [x_optim(4); x_optim(5); x_optim(6)];
R_final=[R_optim T_optim];

end

function [ErrorReproj ]= ObjFunReprojErrLine( x, imgPts, lines,  K)

[err] = calculate_err(x, imgPts, lines,K);
w0 = MM_estimation(err,0,1);
w = MM_estimation(err,w0,2);
ErrorReproj = w .* err' ;


end


function [err] = calculate_err(x, pt, lines,K)

mx =6;
RT = eye(3,4);
RT(1:3,1:3) = rodrigues([x(1) x(2) x(3)]);
RT(1:3,4) = [x(4); x(5); x(6)];
pp= K*RT*lines;
pp= [pp(1,:)./pp(3,:); pp(2,:)./pp(3,:);pp(3,:)./pp(3,:)];
line_kb=lineparasolve(pp)';
n = size(pt, 2);
for i=1:n
    nn=pt(4,i);
    normL = sqrt(line_kb(nn,1)^2 + line_kb(nn,2)^2);
    err(i) = abs(line_kb(nn,1) * pt(1,i)+ line_kb(nn,2) * pt(2,i) + line_kb(nn ,3))/ normL;
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