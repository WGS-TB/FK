
cnf = [ 1 0 1 1 0 0; 0 0 1 1 1 0];
dnf = [ 1 1 1 0 0 0; 0 1 0 1 1 0; 0 0 1 0 0 1];
Check_Conditions(cnf, dnf) % 1 3
Multiple_Check_Conditions(cnf, dnf) 

cnf = [ 0 0 0 0 1; 1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 0 1 1 0 0; 0 1 0 1 0; 0 0 1 1 0];
dnf = [ 1 1 1 0 1];
Check_Conditions(cnf, dnf) % 2 3 4 5
Multiple_Check_Conditions(cnf, dnf) % 2 3 4 5; 1 3 4 5; 1 2 4 5

cnf = [ 0 1 0 1 0 1 0; 0 0 1 1 0 1 0; 0 1 0 1 0 0 1; 0 0 1 1 0 0 1; 1 1 0 0 1 1 0];
dnf = [ 0 1 1 0 0 0 0; 0 0 0 1 1 0 0; 0 0 0 0 0 1 1];
Check_Conditions(cnf, dnf) % 1 3 4 7
Multiple_Check_Conditions(cnf, dnf) 


cnf = [ 1 1 1 0 0 0 0 0; 0 0 0 1 1 1 0 0; 0 0 0 0 0 0 1 1];
dnf = [ 1 1 0 1 0 0 0 0; 1 0 1 0 1 0 0 0];
Check_Conditions(cnf, dnf) % 1 3 4 7
Multiple_Check_Conditions(cnf, dnf) 


cnf = [ 1 1 0; 1 0 1];
dnf = [ 1 0 0];
Check_Conditions(cnf, dnf) % 2 3
Multiple_Check_Conditions(cnf, dnf)

cnf = [ 1 1 0; 1 0 1];
dnf = [ 0 0 1];
Check_Conditions(cnf, dnf) % 3
Multiple_Check_Conditions(cnf, dnf) 


cnf = [ 1 0];
dnf = [ 1 1];
Check_Conditions(cnf, dnf) % 1
Multiple_Check_Conditions(cnf, dnf) 

cnf = [ 1 1 0 0 0; 0 1 1 1 0; 1 0 1 0 0];
dnf = [ 1 0 1 0 0; 1 0 0 1 1; 0 1 1 0 1];
Check_Conditions(cnf, dnf) % 1
Multiple_Check_Conditions(cnf, dnf) 

cnf = [0,1,1,0,0,0;0,1,0,1,0,0];
dnf = [0,0,0,0,1,0];
Check_Conditions(cnf, dnf) % 1
Multiple_Check_Conditions(cnf, dnf) 