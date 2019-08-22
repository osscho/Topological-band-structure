clear all
sigmax = [0 1 ; 1 0];
sigmay = [0 -i ; i 0];
sigmaz = [1 0 ; 0 -1];

A = 1;B = -2; M =2; 
%% gapless (0,pi) at M/2B=-2;
%% gapless (pi,pi) at M/2B=-4;

sampling=20;

for ii = -sampling:sampling
for jj = -sampling:sampling
   
    kx = ii*pi/sampling;
    ky = jj*pi/sampling;
      
    hab = -2*(cos(kx)*exp(-i*pi/4)+cos(ky)*exp(i*pi/4));
    hab_star = -2*(cos(kx)*exp(i*pi/4)+cos(ky)*exp(-i*pi/4));
   
    %% Hamiltonian
    Hamilt = [0, hab; hab_star,0];
    
    [Evectors,Evalues] = eig(Hamilt);

    Eval_mat_hig(ii+sampling+1,jj+sampling+1) = Evalues(1,1);
    Eval_mat_low(ii+sampling+1,jj+sampling+1) = Evalues(2,2);

    
    Evec_mat_left = Evectors(:,1);
    Evec_mat_right= Evectors(:,2);
   

end
end


figure
subplot(221)
hold
shading faceted
surf(Eval_mat_hig,'edgecolor','none')
surf(Eval_mat_low,'edgecolor','none')

light
camlight('right')
camproj('orthographic')


