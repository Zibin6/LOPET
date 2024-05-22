function wi = MM_estimation(err,w0,initial)

if initial==1
err_median = median(err) ;
MAD = median(abs(err-err_median));
err = err/(1.4826.*MAD);

for i=1:size(err,2)
    if abs(err(i))<=4.685
        wi(i)=(1-(err(i)/4.685).^2).^2;
    else
        wi(i)=0;
    end
end

wi = wi';

else

sigma_S = sqrt(mean(w0'.*err.*err)/0.199);
err = err/sigma_S;


for i=1:size(err,2)
    if abs(err(i))<=4.685
        wi(i)=(1-(err(i)/4.685).^2).^2;
    else
        wi(i)=0;
    end
end


wi = wi'/sum(wi);
end


end