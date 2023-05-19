function[W] = coupling_matrix (Q, n, m)
 k = 0;
   for i = 1:m
   l = 0;
     for j = 1:m
     W(i,j) = sum(sum(Q(k+1 : k+n(i), l+1 : l+n(j))))/n(i);
     l = l + n(j);
     end
   k = k + n(i);
 end

 end
