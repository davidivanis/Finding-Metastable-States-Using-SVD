%stohasti?ka matrica s tri metastabilna

A=[0.35 0.15 0.5 0 0 0 0 0 0 ;
   0.2 0.75 0.05 0 0 0 0 0 0 ;
   0.3 0.5 0.2 0 0 0 0 0 0 ;
   0 0 0 0.24 0.57 0.19 0 0 0;
   0 0 0 0.32 0.36 0.32 0 0 0;
   0 0 0 0.2 0.7 0.1 0 0 0;
   0 0 0 0 0 0 0.2 0.6 0.2;
   0 0 0 0 0 0 0.2 0.55 0.25;
   0 0 0 0 0 0 0.1 0.85 0.05];
   
%matrica perturbacije, na dijagonalnim blokovima uniformno na (-1,1) okolo na (0,1)
E=rand(9,9);
for i=1:3
  E((i-1)*3+1:i*3,(i-1)*3+1:i*3)=-1+2*rand(3,3);
end

%perturbiramo A sa E i vratimo ju na stohasti?ku
eps=0.5*1e-1;
M=A+eps*E;

for i=1:9
  M(i,:)=M(i,:)/sum(M(i,:));
end

%svojstveni vektori i vrijednosti od M, prikazani grafi?ki

[V,S,W]=eig(M);
V=(1/V(1,1))*V; %da mi u prvom svojstvenom budu sve jedinice, zbog prikaza
x=[1:9];
xi=[x(1):0.001:x(end)]; %da plotam svojstvene kao glatki graf a ne to?ke
for i=1:3
  figure
  plot(x,V(:,i),'b*');
  vid=interp1(x,V(:,i),xi,'spline');
  hold on
  plot(xi,vid,'r','linewidth',1.2);
  line([1,9],[0,0],'color','b');
  line([3.5,3.5],[-2,2],'color','b','linestyle','-.');
  line([6.5,6.5],[-2,2],'color','b','linestyle','-.');
  if i==1
   plot(2.25,0.3,'+k','linewidth',1.2);
   plot(5,0.3,'+k','linewidth',1.2);
   plot(7.75,0.3,'+k','linewidth',1.2);
  
  else 
    if V(1,i)<0
      line([2.15,2.35],[0.3,0.3],'color','k','linewidth',1.2);
    else
      plot(2.25,0.3,'+k','linewidth',1.2);
    end
    if V(4,i)<0
      line([4.9,5.1],[0.3,0.3],'color','k','linewidth',1.2);
    else
      plot(5,0.3,'+k','linewidth',1.2);
    end
    if V(7,i)<0
      line([7.65,7.85],[0.3,0.3],'color','k','linewidth',1.2);
    else
      plot(7.75,0.3,'+k','linewidth',1.2);  
    end
  end
  
  axis([1,9,-2,2]);
  ime=strcat("sv_vek_",mat2str(i));
  saveas(i,ime,'jpg');
end



%ista stvar ali sa singularnim vektorima


[U,D,UU]=svd(M);
for i=1:3
  figure
  plot(x,U(:,i),'b*');
  vid=interp1(x,U(:,i),xi,'spline');
  hold on
  plot(xi,vid,'r','linewidth',1.2);
  line([1,9],[0,0],'color','b');
  line([3.5,3.5],[-2,2],'color','b','linestyle','-.');
  line([6.5,6.5],[-2,2],'color','b','linestyle','-.');
  
  if U(1,i)<0
      line([2.15,2.35],[0.3,0.3],'color','k','linewidth',1.2);
  else
      plot(2.25,0.3,'+k','linewidth',1.2);
  end
  if U(4,i)<0
      line([4.9,5.1],[0.3,0.3],'color','k','linewidth',1.2);
  else
      plot(5,0.3,'+k','linewidth',1.2);
  end
  if U(7,i)<0
      line([7.65,7.85],[0.3,0.3],'color','k','linewidth',1.2);
  else
      plot(7.75,0.3,'+k','linewidth',1.2);  
  end
  
  axis([1,9,-2,2]);
  ime=strcat("sing_vek_",mat2str(i));
  saveas(i+3,ime,'jpg');
end


%ispermutiram A, i na njoj testiram kod. Prikazujem po?etnu A, njenu permutaciju i rezultat algoritma grafi?kill

figure
imagesc(M);
colormap(flipud(hot));
colorbar;
ime="M_prije";
saveas(7,ime,'jpg');

Perm=eye(9);
Perm=Perm(randperm(9),:);
M_permutirana=Perm*M*Perm';

figure
imagesc(M_permutirana);
colormap(flipud(hot));
colorbar;
ime="M_permutirana";
saveas(8,ime,'jpg');

[m,n,P]=SVD_metastable(M_permutirana,9,1,9,0.7,0,0,eye(9));

M_poslije=P*M_permutirana*P';

figure
imagesc(M_poslije);
colormap(flipud(hot));
colorbar;
ime="M_poslije";
saveas(9,ime,'jpg');

figure
e=eig(M);
s=sort(abs(e),'descend');
plot(1:9,s(1:9),'.b','MarkerSize',10);
axis([0,10,0,1]);
ime="sv_vr_M";
saveas(10,ime,'jpg');

figure
e=diag(D);
s=sort(abs(e),'descend');
plot(1:9,s(1:9),'.b','MarkerSize',10);
axis([0,10,s(9)-0.1,s(1)+0.1]);
ime="sing_vr_M";
saveas(11,ime,'jpg');

C=coupling_matrix(M_poslije,n,m);