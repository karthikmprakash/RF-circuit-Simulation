%****************************************
% ELG 7132D. Assignment 2
% Student name: Karthik Mysore Prakash
% Student ID Numeber: 300037618
% Date: October 20 2017
%****************************************



%% The first step is to initialize the global data structures and add some default paths to the Matlab engine
clc
close all
clear all 
xDC=0; 
init_globals

%% Defining the file describing the circuit netlist
circuit_file_name = './netlists/Assign-2-common-emitter.sp';

%% Obtaining the Initial Guess x(0)
my_circuit_struct = mna_parse_circuit(circuit_file_name);
X_0_init_guess = my_circuit_struct.xic
xi= X_0_init_guess;

%% Obtaining the MNA 'G' Matrix
G=full(my_circuit_struct.G)

C = full(my_circuit_struct.C)

%% Comuting the DC Source Vector 
bDC = mna_compute_source_vector(my_circuit_struct,0)

%% Computing the non linear function vector f(x(t)) and Jacobian matrix of the non linear part pdf(x)/pdx 
i=0;
[fx b Jf] = mna_compute_circuit(my_circuit_struct,0,xi)
phi_xi=(G*xi)+ fx -bDC
e=power(10,-14);
phi_xi_norm=norm(phi_xi);
while phi_xi_norm>e
Jxi=(G+Jf);
inverse_Jxi= inv(Jxi);
xi_new=xi-(inverse_Jxi*phi_xi);
xi= xi_new;
[fx b Jf] = mna_compute_circuit(my_circuit_struct,0,xi)
phi_xi=(G*xi)+fx-bDC
phi_xi_norm=norm(phi_xi);
i=i+1;
xDC=xi
end

b=(C*Jf)+(G*xi)+fx;


