clc
clear all
sigmax = [0 1 ; 1 0];
sigmay = [0 -i ; i 0];
sigmaz = [1 0 ; 0 -1];

t2=0.2;
sampling=100;
layers=20;

figure
hold on

for ii = -sampling:sampling

              
    kx = ii*pi/sampling;

        
   %%
     
    Hnn1=[0, 2*cos(kx)*(1+i);2*cos(kx)*(1-i),0];
      %%
       I = eye(layers);
       Hnnc1=kron(I,Hnn1);
      %%
%       Hnn2=exp(-i*ky*(jj-jjj-1))*[0,1-i;1+i,0];
     
      A23=[0,1+i];
      B23=repmat(A23,[1,layers]);
      Hnnup2=diag(B23,1);
      Hnnup22=Hnnup2(1:2*layers,1:2*layers);
      
      A14 = [1-i,0];
      B14=repmat(A14,[1,layers-1]);
      Hnnup3=diag(B14,3);
      Hnnup33=Hnnup3(1:2*layers,1:2*layers);
           
      
      %%
%       Hnn3=exp(-i*ky*(jj-jjj+1))*[0,1-i;1+i,0];
      A32=[0,1-i];
      B32=repmat(A32,[1,layers]);
      Hnndown2=diag(B32,-1);
      Hnndown22=Hnndown2(1:2*layers,1:2*layers);
      
      A41=[1+i,0];
      B41=repmat(A41,[1,layers-1]);
      Hnndown3=diag(B41,-3);
      Hnndown33=Hnndown3(1:2*layers,1:2*layers);
      
     Hnn=Hnnc1+Hnnup22+Hnnup33+Hnndown22+Hnndown33;
     %% H next nearest neighbor
%       Hnnn1=i*2*sin(kx)*exp(-i*ky*(jj-jjj-1))*[-1,0;0,1];
%       Hnnn2=i*2*sin(kx)*exp(-i*ky*(jj-jjj+1))*[1,0;0,-1];

      C3=[i*sin(kx),-i*sin(kx)];
      C33=repmat(C3,[1,layers]);
      Hnnnup=diag(C33,2);
      Hnnnup1=Hnnnup(1:2*layers,1:2*layers);
      
      %%
      
     D3=[-i*sin(kx),i*sin(kx)];
     D33=repmat(D3,[1,layers]);
      Hnnndown=diag(D33,-2);
      Hnnndown1=Hnnndown(1:2*layers,1:2*layers);
      Hnnn=Hnnnup1+Hnnndown1;
      %%
     Hamilt=Hnn+t2*Hnnn;
 
%         
     
   [Evectors,Evalues] = eig(Hamilt);
   

 plot(kx,diag(Evalues),'.')

 



end

 
