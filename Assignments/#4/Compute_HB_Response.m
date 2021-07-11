function [frequency,Xo,Xc,Xs,THD] = Compute_HB_Response(nTimePoints,node_name)
%****************************************
% ELG 7132D. Assignment 4
% Student name: Karthik Mysore Prakash
% Student ID Numeber: 300037618
% Date: Nov 18 2017
%Function for computing transient response using HB
%****************************************
init_globals
circuit_file_name = './netlists/Assign-2-common-emitter.sp';
  my_circuit_struct = mna_parse_circuit('./netlists/Assign-2-common-emitter.sp');
  
  node_value = mna_get_var_index(my_circuit_struct.circuit_id,node_name,'node');
  nv=node_value;
  H=nTimePoints;
  N= length(my_circuit_struct.xic); 
  f=5e6;
  T=1/f;
  w0=2*pi*f;
  
  %Checking for odd or even
 if (rem (nTimePoints,2) == 0)
    k = (nTimePoints-2)/2;
    k=k+1;
    O=zeros(N,((2*k)-1)*N); % for Y_bar
else
    k = (nTimePoints-1)/2;
    O=zeros(N,2*k*N);
    
end
  G=full(my_circuit_struct.G);
  C=full(my_circuit_struct.C);
  
 
  %% ****************** Creating B(T) *************************************
   DeltaT = T/H;
TimePoints = [0:DeltaT:T-DeltaT]';
bt = zeros(N,length(TimePoints));
for t=1:length(TimePoints)
    bt(:,t) = mna_compute_source_vector(my_circuit_struct, TimePoints(t));
end
     bt_transpose = bt';
bt_bar = bt_transpose(:);
%Creating Gamma_bar
  gamma_inv=makeRealGammaInv(nTimePoints);
  gamma= inv(gamma_inv);
  gamma_bar=kron(eye(N),gamma);
  gamma_bar_inv=kron(eye(N),gamma_inv);
  % Creating P
  eP=zeros(N,1);
  ePT=zeros(1,N);
  IH=eye(H);
  P=0;
  for i=1:1:N
      eP(i,1)=1;
      ePT=eP';
      P1=kron(eP,IH);
      P2=kron(P1,ePT);
      P=P+P2;
      eP(i,1)=0;
  end
  PT=P';
    BT_bar=PT*gamma_bar*bt_bar;
  BDc_bar=zeros(N*H,1);
   BDc_bar(1:N,1)=BT_bar(1:N,1);
   
   BAc_bar=zeros(N*H,1);
   BAc_bar(N:end,1)=BT_bar(N:end,1);

  %% *************** Generating initial X_bar *****************************
  % Initial guess
  X_0_init_guess = my_circuit_struct.xic;
  xi=X_0_init_guess;
  G=full(my_circuit_struct.G);
bDC = mna_compute_source_vector(my_circuit_struct,0); %deriving the DC source value from the netlist
i=0;
[fx b Jf] = mna_compute_circuit(my_circuit_struct,0,xi); % deriving fx and Jf from the netlist
phi_xi=(G*xi)+ fx-bDC; %error of the guess
e=10e-14;
phi_xi_norm=norm(phi_xi); % norm of the error
while phi_xi_norm>=e, % to check for minimum error
Jxi=(G+Jf); % Jacobian matrix formation
inverse_Jxi= inv(Jxi);
xi_new=xi-(inverse_Jxi*phi_xi);% new guess is found depending on the equation
xi= xi_new; % initial guess is being replaced by the new guess
[fx b Jf] = mna_compute_circuit(my_circuit_struct,0,xi); % fx and Jf is found for the new guess.
phi_xi=(G*xi)+fx-bDC;  %error for the new guess
phi_xi_norm=norm(phi_xi); %norm of the phi for new guess is found. If phi satisfy the given 'while' conditon the loop will be running until the value of phi is less than 'e'
 end
xDC=xi; %the value of xi for which it comes out of while loop is our desired xDC
% end of initial guess
  
% X_bar
X_bar=zeros(N*H,1);

X_bar(1:N,1)=xDC(1:N,1);

% End of X_Bar


  %% **************** Compute HB response Start ***************************
  %Creating Y bar
  I2=eye(2);
  T2=[0 1;-1 0];
  OT=O';
  e=zeros(k,1);
  eT=e';
  %Creating the sum
  YB_sum=0;
  for m=1:1:k
      e(m,1)=1;
      eT=e';
      L=e*eT;
      L1=kron(I2,G);
      mwT=(m*w0).*T2;
      L2=kron(mwT,C);
      L3=L1+L2;
      L4=kron(L,L3);
      YB_sum=YB_sum+L4;
      e(m,1)=0;
  end
 %generating YB_sum for even or odd timepoints
   if (rem (nTimePoints,2) == 0)
    k = (nTimePoints-2)/2;
    YB_Sum=YB_sum(1:N*(H-1),1:N*(H-1));
else
    k = (nTimePoints-1)/2;
    YB_Sum=YB_sum;
end
  
  Y_bar=[G O;OT YB_Sum];
  
 %% Creating F_bar(X_bar)
  % Creating P
  eP=zeros(N,1);
  ePT=zeros(1,N);
  IH=eye(H);
  P=0;
  for i=1:1:N
     eP(i,1)=1;
     ePT=eP';
      P1=kron(eP,IH);
      P2=kron(P1,ePT);
      P=P+P2;
      eP(i,1)=0;
  end
  PT=P';
  %Creating Gamma_bar
  gamma_inv=makeRealGammaInv(nTimePoints);
  gamma= inv(gamma_inv);
  gamma_bar=kron(eye(N),gamma);
  gamma_bar_inv=kron(eye(N),gamma_inv);
 
  
  c=0;
  for a=0.01:0.01:1 %to implement the split BT_bar calculation
  e=10e-14;
  phinm=1;
  while phinm >e,
  %generating the sum inside F_bar(X_bar)
   FX1=gamma_bar_inv*P*X_bar;
   IN=eye(N);
   ef=zeros(H,1);
   FX2=0;
   for j=1:1:H
       ef(j,1)=1;
       efT=ef';
       f1=kron(IN,efT);
       f2=f1*FX1;
       [fx b Jf] = mna_compute_circuit(my_circuit_struct,0,f2);% deriving fx and Jf from the netlist
  FX=kron(fx,ef);
  FX2=FX2+FX;
  ef(j,1)=0;  
   end
   FX_bar=PT*gamma_bar*FX2;
   
   %% Creating Jacobian
    ej=zeros(H,1);
      jnode=0;
      for j=1:1:H
         ej(j,1)=1;
   ejt=ej';
   xj=(kron(IN,ejt)*gamma_bar_inv*P*X_bar);
   [fx b Jf] = mna_compute_circuit(my_circuit_struct,0,xj);
   j1=gamma*ej*ejt*gamma_inv; % j1 becomes ill conditioned as Lm,k is ill conditioned
   j2=kron(Jf,j1);
   jnode=jnode+j2;
   ej(j,1)=0;   
      end
       
   del_FX_bar=PT*jnode*P;
   J=Y_bar+del_FX_bar; %Jacobian Matrix
   Jinv=inv(J);        %May Contain some portion of error being ill conditioned
          phiX_bar=(Y_bar*X_bar)+FX_bar-BDc_bar-a*BAc_bar;
          X1_bar=X_bar-(Jinv*phiX_bar); %updating the value of X_bar
          X_bar=X1_bar;
              phinm=norm(phiX_bar,2);
              phi=phinm;
  end
  
  end
   % Generating Frequency Value  
frequency=0;
  
   for i=0:1:k
   frequency(i+1,1)=i*5e6;
   end
   
   %Generating frequency components 
      Xo=X_bar(nv,1);
      Xc(:,1)=X_bar((nv+N):2*N:N*H-N,1);
      Xs(:,1)=X_bar((nv+(2*N)):2*N:N*H,1);
      
   %% Computing total harmonic distortion
   THD=0;
   sum=0;
   fComp=0;
   for i=2:1:k
       sum1=sqrt((Xc(i,1)^2)+(Xs(i,1)^2));
      sum=sum+sum1;
      fComp(i,1)=sum1;
   end
   den=sqrt((Xc(1).^2+ Xs(1).^2));
   THD=log10(sum1/den)
       

end
function GammaInv = makeRealGammaInv(nTimePoints)
if (rem (nTimePoints,2) == 0)
    % even number of points
    % truncate one positive frequency component
    K = (nTimePoints-2)/2;

    BC = [1 0]; BS = [0 1]; B1 = [1:K];
    M = [0:nTimePoints-1]'*[B1]; 
    Ma = [0:nTimePoints-1]'*(K+1);  
    MC = cos(M*2*pi/(nTimePoints));  
    MS = sin(M*2*pi/(nTimePoints));

    MCa = [cos(Ma*2*pi/(nTimePoints))];
    MSa = zeros(nTimePoints,1); 
    MC1 = [kron(MC,BC) MCa];
    MS1 = [kron(MS,BS) MSa];
    GammaInv = [ones(nTimePoints,1) MC1+MS1];
else    % odd number of points
   
    K = (nTimePoints -1)/2;
    BC = [1 0];
    BS = [0 1];
    B1 = [1:K];
    B2C = kron(B1,BC);
    B2S = kron(B1,BS);
    M = [0:nTimePoints-1]'*[B1];

    MC = cos(M*2*pi/(nTimePoints)); MS = sin(M*2*pi/(nTimePoints));
    MC1 = kron(MC,BC);  MS1 = kron(MS,BS);
    GammaInv = [ones(nTimePoints,1) MC1+MS1];

end
end