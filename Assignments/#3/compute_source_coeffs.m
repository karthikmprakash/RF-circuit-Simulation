  function [B0,Bc,Bs] = compute_source_coeffs(n_Time_Points)
 %****************************************
% ELG 7132D. Assignment 3
% Student name: Karthik Mysore Prakash
% Student ID Numeber: 300037618
% Date: November 02 2017
%Function for computing Source Coefficients
% It is needed to clear the workspace before using this function
%****************************************

gamma_inv=zeros(n_Time_Points); % Initializing the gamma inverse matrix with zeros and size entered
H=length(gamma_inv);



init_globals
%  Defining the file describing the circuit netlist
circuit_file_name = './netlists/tuned_amplifier_1.sp';
my_circuit_struct = mna_parse_circuit(circuit_file_name);


T = 5e-9;          %The timeperiod of the input
Time = linspace(0,T-T/H,H);

R=1;    % Variable assigned to Rows
i=1;    % Variable that keeps a count of the time iterations
        % x(H-1,1)=0;

%%
if(mod(H,2) == 0) %To check for the odd or even time points
    for t=0:T/H:(T-(T/H)) %For even time points
            b=mna_compute_source_vector(my_circuit_struct,t); %The index of the time sample
            gamma_inv(R,1)= 1; %The First column of the Gamma INverse matrix should be one
            k=1; % Variable in the Fourier series
    
                 for C=2:2:H-2
                     gamma_inv(R,C)= cos(k*i*2*pi/H);
                     C=C+1; %Column Increment
                     gamma_inv(R,C)=sin(k*i*2*pi/H);
                     k=k+1;
                 end
       C=C+1;
       gamma_inv(R,C)=cos(k*i*2*pi*(1/H)); % Computing the cos for the odd or the H-1 term
       time(R,1)=t; %Keeping track of time for plotting 
       i=i+1;
      
      x(R,1)=b(15,1); % Computing the input source vector for the time domain tracing
       R=R+1; % Row increment
    end
    
else 
        for t=0:T/H:(T-(T/H)) 
            b=mna_compute_source_vector(my_circuit_struct,t);                %The index of the time sample
            gamma_inv(R,1)= 1;
            k=1;
    
                 for C=2:2:H-1
                     gamma_inv(R,C)= cos(k*i*2*pi/H);
                     C=C+1; %Column Increment
                     gamma_inv(R,C)=sin(k*i*2*pi/H);
                     k=k+1;
                 end
      
             time(R,1)=t; %Keeping track of time for plotting
            i=i+1;
       
            x(R,1)=b(15,1);% Computing the input source vector for the time domain tracing
            R=R+1;% Row increment
        end     
end


%% Computing the output and tracing

gamma_trsp=gamma_inv';
gamma=(2/H).*gamma_trsp;
gamma_ROW1=gamma(1,1:end).*(1/2);
gamma(1,:)=gamma_ROW1;
X=gamma*x;

%The coefficients of the source vector 
B0=X(1,1);
Bs=X(3:2:length(x));
Bc=X(2:2:length(x));

x_new=gamma_inv*X; %using X(T) to find the initial sample

%Tracing the Outputs obtained from Gamma Inverse matrix 
%Computation and the time domain computation
figure(1)
subplot(2,1,1)
plot(time,x,'LineWidth',4)
title('Output from the netlist');
xlabel('Time in secs');
ylabel('Aplitude');
set(gca,'fontsize',20)
%  hold on
subplot(2,1,2)
plot(time,x_new,'--g','LineWidth',4)
title('Time Domain Output');
xlabel('Time in secs');
ylabel('Aplitude');
set(gca,'fontsize',20)
grid on 

figure(2)
plot(time,x,'LineWidth',4)
hold on
plot(time,x_new,'--g','LineWidth',4)
grid on
title('Super Imposed Trace');
xlabel('Time in secs');
ylabel('Aplitude');
set(gca,'fontsize',20)
legend('Netlist Output','Computed output')

mna_destroy_circuit(my_circuit_struct)
 end