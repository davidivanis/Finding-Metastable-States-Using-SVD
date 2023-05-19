function[m, n, P] = SVD_metastable (B, nB , a, b, r, m, n, P)
 %izracun prva dva svojstvena vektora
 [U,S,V] = svds(B ,2);
 %sortiranje drugog singularnog vektora i kreiranje P
 [u2 , I] = sort(U(: ,2));
 Pp = eye(size(u2 ,1));
 Pp = Pp(I ,:);
 %primjena permutacije na matricu B
 B2 = Pp*B*(Pp)';
 %identifikacija dva potencijalna bloka
 ind = find(diff(sign(u2)));

 D = eye(nB);
 D(a:b,a:b)=Pp;
 P = D*P;
 %izdvojimo dobivene potencijalne blokove
 B_prvi = B2 (1:ind ,1: ind);
 B_drugi= B2(ind +1:size(B ,1) , ind +1:size(B ,1));
 %racunanje 1-norme blokova
 normaB1 = sum(sum(B_prvi))/size(B_prvi ,1);
 normaB2 = sum(sum(B_drugi))/size(B_drugi ,1);
 if normaB1 > r & normaB2 > r
 [m, n, P] = SVD_metastable (B_prvi , nB , a, a+ind -1, r, m, n, P);
 [m, n, P] = SVD_metastable (B_drugi , nB , ind+a, b, r, m, n, P);
 return

 else
 m = m + 1;
 n(m) = size(B ,1);
 return;
 end
 end
