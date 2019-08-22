clear all
format long
sigmax = [0 1 ; 1 0];
sigmay = [0 -i ; i 0];
sigmaz = [1 0 ; 0 -1];

A = 1;B =-2; M =1; 

s=200;

Hamilt=zeros(2,2,s);
dH_dkx=zeros(2,2,s);
dH_dky=zeros(2,2,s);
Curvature= zeros(s,s);

kx=linspace(-pi,pi,s);
ky=linspace(-pi,pi,s);
 
for ii = 1:length(kx)
    for jj= 1:length(ky)
   
	dx = A * sin(kx(ii));
    dy = A * sin(ky(jj));
    dz = M + 2*B*( 2-cos(kx(ii))-cos(ky(jj)) );
      %% Hamiltonian
    Hamilt(:,:,jj) = dx*sigmax + dy*sigmay + dz*sigmaz; 
   
    [Evectors,Evalues] = eig(Hamilt(:,:,jj));
     
     E_1(ii,jj) = Evalues(1,1);
     E_2(ii,jj) = Evalues(2,2);
%    
    dH_dkx(:,:,ii)=[2*B*sin(kx(ii)), A*cos(kx(ii)); A*cos(kx(ii)),  -2*B*sin(kx(ii))];      
    dH_dky(:,:,jj)=[2*B*sin(ky(jj)), -i*A*cos(ky(jj)); i*A*cos(ky(jj)), -2*B*sin(ky(jj))];
% 
    Num1= Evectors(:,1)'*dH_dkx(:,:,ii)*Evectors(:,2);
    Num2 = Evectors(:,2)'*dH_dky(:,:,jj)*Evectors(:,1);

    Num3= Evectors(:,1)'*dH_dky(:,:,jj)*Evectors(:,2);
    Num4= Evectors(:,2)'*dH_dkx(:,:,ii)*Evectors(:,1);

    Numerator= Num1*Num2-Num3*Num4;
    dellamda=(Evalues(1,1)-Evalues(2,2))^2;
     
    Curvature(ii,jj)=i*Numerator./dellamda/(s-1)/(s-1)*pi*pi;
    Chern = sum(sum(Curvature))/4/4;
   end
end

Chern

subplot(221)
hold on
surf(E_1);
surf(E_2);
hold off

subplot(222)
surf(real(Curvature))

subplot(223)
hold on
surf(real(Hamilt(:,:)))
surf(imag(Hamilt(:,:)))
hold off
subplot(224)
hold on
surf(dH_dkx(:,:))
surf(imag(dH_dky(:,:)))
hold off

