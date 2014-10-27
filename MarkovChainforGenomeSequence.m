clear all;
close all;
clc;
% Markov Chain for genome sequence
%     A     C     G     T
P=[0.15 0.166 0.1875 0.2; ... %Probability transition matrix
    0.35 0.334 0.3125 0.3; ...
    0.35 0.3125 0.3125 0.3; ...
    0.15 0.1875 0.1875 0.2]

% Probability of the following three sequences 
%P(AAGCCCTGGCAATTCAG)=
P1=P(1,1)*P(3,1)*P(2,3)*P(2,2)*P(2,2)*P(4,2)*P(3,4)*P(3,3)*P(2,3) ...
    *P(1,2)*P(1,1)*P(4,1)*P(4,4)*P(2,4)*P(1,2)*P(3,1)
%P(ACCTGCGCCGTATATTA)=
P2=P(2,1)*P(2,2)*P(4,2)*P(3,4)*P(2,3)*P(3,2)*P(2,3)*P(2,2)*P(3,2) ...
    *P(4,3)*P(1,4)*P(4,1)*P(1,4)*P(4,1)*P(4,4)*P(1,4)
%P(GGCTCTCCAAGCCTTAT)=
P3=P(3,3)*P(2,3)*P(4,2)*P(2,4)*P(4,2)*P(2,4)*P(2,2)*P(1,2)*P(1,1) ...
    *P(3,1)*P(2,3)*P(2,2)*P(4,2)*P(4,4)*P(1,4)*P(4,1)
%% 
[V,D]=eig(P); %eig returns eigenvectors in their respective columns and their associated eigenvalues
Aproportion=V(1,1)/(sum(V(:,1))); %calculates the probability of encountering A in the genome
Cproportion=V(2,1)/(sum(V(:,1))); %calculates the probability of encountering C in the genome
Gproportion=V(3,1)/(sum(V(:,1))); %calculates the probability of encountering G in the genome
Tproportion=V(4,1)/(sum(V(:,1))); %calculates the probability of encountering T in the genome

% Calculate the mean recurrence time for "A" in a genome sequence
%Setup the transition probability matrix
%     A     C      G    T
P = [0.15 0.166 0.1875 0.2; ... %probability transition 	matrix
    0.35 0.334 0.3125 0.3; ...
    0.35 0.3125 0.3125 0.3; ...
    0.15 0.1875 0.1875 0.2];
%Extract the submatrix of transient state transition 	probabilities
T=P([2:end],[2:end]);
% Calculate the fundamental matrix, W
W=inv(eye(3)-T);
% Output the vector of mean time to absorption by summing along the rows
tau=sum(W,1);
% Calculating the mean number of letters between two A's
% P_AA=P(1,1);
% P_AC=P(2,1);
% P_AG=P(3,1);
% P_AT=P(4,1);
% tau_C=tau(1);
% tau_G=tau(2);
% tau_T=tau(3);
% tau_A=P_AA*1+P_AC*(1+tau_C)+P_AG*(1+tau_G)+P_AT*(1+tau_T)
tau_A=P(1,1)*1+P(2,1)*(1+tau(1))+P(3,1)*(1+tau(2))+P(4,1)*(1+tau(3));	
