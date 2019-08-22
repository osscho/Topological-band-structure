%% The aim is to simulate Berry curvature according to topological band Theory (Haldane model)
%% Output: Berry curvation = 1 
clear all
sigmax = [0 1 ; 1 0];
sigmay = [0 -i ; i 0];
sigmaz = [1 0 ; 0 -1];

A = 1;B = -2; M =1; 
%% gapless (0,pi) at M/2B=-2;
%% gapless (pi,pi) at M/2B=-4

sampling =100;

for ii = -sampling:sampling
for jj = -sampling:sampling
   
    kx = ii*pi/sampling;
    ky = jj*pi/sampling;
      
    dx = A * sin(kx);
    dy = A * sin(ky);
    dz = M + 2*B*( 2-cos(kx)-cos(ky) );

    %% Hamiltonian
    Hamilt = dx*sigmax + dy*sigmay + dz*sigmaz;
    
    dkx=[2*B*sin(kx), A*cos(kx); A*cos(kx),  -2*B*sin(kx)];     
    dky=[2*B*sin(ky), -i*A*cos(ky); i*A*cos(ky), -2*B*sin(ky)];
    
    [Evectors,Evalues] = eig(Hamilt);

    Eval_mat_hig(ii+sampling+1,jj+sampling+1) = Evalues(1,1);
    Eval_mat_low(ii+sampling+1,jj+sampling+1) = Evalues(2,2);

    
    Evec_mat_left = Evectors(:,1);
    Evec_mat_right= Evectors(:,2);

    Num1 = Evec_mat_left'*dkx*Evec_mat_right;
    Num2 = Evec_mat_right'*dky*Evec_mat_left;
   
    Num3= Evec_mat_left'*dky*Evec_mat_right;
    Num4= Evec_mat_right'*dkx*Evec_mat_left;

    Numerator= Num1*Num2-Num3*Num4;
    Denominator=(Evalues(1,1)-Evalues(2,2)).^2;
    
   Curvature(ii+sampling+1,jj+sampling+1)= i*Numerator./Denominator/(2*sampling)/(2*sampling)*pi*pi*pi*pi;

  
end
end

Sum=sum(sum(Curvature))/4/4
figure
subplot(221)
hold
shading faceted
surf(Eval_mat_hig,'edgecolor','none')
surf(Eval_mat_low,'edgecolor','none')

light
camlight('right')
camproj('orthographic')


subplot(222)
surf(real(Curvature),'edgecolor','none')
light
camlight('right')
camproj('pers')
grid off



