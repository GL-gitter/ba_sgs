A = load("MA_Jonas\matlab.mat");
B = load("MA_Jonas\matlab2.mat");

MatA = A.A;
MatB = B.A;

MatC = MatA-MatB;


MatC1 = MatA+MatB;