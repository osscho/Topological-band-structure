clear all
clc

sigmax = [0 1 ; 1 0]
sigmay = [0 -i ; i 0]
sigmaz = [1 0 ; 0 -1]

A = 1;
B = -1;
M = 2;

sampling = 20;

for ii = -sampling:sampling
for jj = -sampling:sampling

    kx = ii*pi/sampling;
    ky = jj*pi/sampling;
    dx = A * sin(kx);
    dy = A * sin(ky);
    dz = M + 2*B*( 2-cos(kx)-cos(ky) );

    Hamilt = dx*sigmax + dy*sigmay + dz*sigmaz ;
    
    [Evectors,Evalues] = eig(Hamilt);

    Eval_mat_low(ii+sampling+1,jj+sampling+1) = Evalues(1,1);
    Eval_mat_hig(ii+sampling+1,jj+sampling+1) = Evalues(2,2);
    
end
end

figure
hold

surf(Eval_mat_hig)
surf(Eval_mat_low)

