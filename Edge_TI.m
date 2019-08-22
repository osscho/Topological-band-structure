clc
clear all
sigmax = [0 1 ; 1 0];
sigmay = [0 -i ; i 0];
sigmaz = [1 0 ; 0 -1];

format long

A = 1;B = -2; M =4; 
t2=0.2;
sampling=20;
layers=20;

figure
hold on

for ii = -sampling:sampling
              
    kx = ii*pi/sampling;
        
   %%
     
    Hnn1=[0,sin(kx);sin(kx),0];
      %%
       I = eye(layers);
       Hnnc1=kron(I,Hnn1);
      %%

     
      A23=[0,0.5];
      B23=repmat(A23,[1,layers]);
      Hnnup2=diag(B23,1);
      Hnnup22=Hnnup2(1:2*layers,1:2*layers);
      
      A14 = [-0.5,0];
      B14=repmat(A14,[1,layers-1]);
      Hnnup3=diag(B14,3);
      Hnnup33=Hnnup3(1:2*layers,1:2*layers);
           
      
      %%
      A32=[0,0.5];
      B32=repmat(A32,[1,layers]);
      Hnndown2=diag(B32,-1);
      Hnndown22=Hnndown2(1:2*layers,1:2*layers);
      
      A41=[-0.5,0];
      B41=repmat(A41,[1,layers-1]);
      Hnndown3=diag(B41,-3);
      Hnndown33=Hnndown3(1:2*layers,1:2*layers);
      
     Hnn=Hnnc1+Hnnup22+Hnnup33+Hnndown22+Hnndown33;
    %% H next nearest neighbor
%%
       
Hnnndiag=[M+2*B*(2-cos(kx)),0;0,-M-2*B*(2-cos(kx))];
Hnnn1=kron(I,Hnnndiag);

%%


      C3=[-1,1];
      C33=repmat(C3,[1,layers]);
      Hnnnup=diag(C33,2);
      Hnnnup1=Hnnnup(1:2*layers,1:2*layers);
      
      %%
      
     D3=[-1,1];
     D33=repmat(D3,[1,layers]);
      Hnnndown=diag(D33,-2);
      Hnnndown1=Hnnndown(1:2*layers,1:2*layers);
      Hnnn=Hnnn1+Hnnnup1+Hnnndown1;
      %%
     Hamilt=Hnn+Hnnn;
 %         
     
   [Evectors,Evalues] = eig(Hamilt);
   

 plot(kx,diag(Evalues),'.','MarkerEdgeColor','b') 


end




