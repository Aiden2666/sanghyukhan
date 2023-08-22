function tm = tmean (x, pct)
% x is the array to avarage
% pct is the percentage of values that we want to keep (0-99)

l      = length(x);
x_sort = sort(x);

lo_ind = int8(floor ( ((100-pct) /2 ) *l/100 ) )
hi_ind = int8(l-lo_ind)

x_sort(hi_ind:end)=[];
x_sort(1:lo_ind)=[]

tm = mean(x_sort);



