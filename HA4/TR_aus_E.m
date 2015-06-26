function [T1,R1,T2,R2]=TR_aus_E(E)
% In dieser Funktion sollen die moegleichen euklidischen Transformationen
% aus der essentiellen Matrix extrahiert werden

% Calculate SVD of E
[U,S,V]=svd(E);

R_zp=[0,-1,0;
      1,0,0;
      0,0,1];
R_zm=[0,1,0;
      -1,0,0;
      0,0,1];
  
R1=U*R_zp'*V';
R2=U*R_zm'*V';
T1=inverse_dachoperator(U*R_zp*S*U');
T2=inverse_dachoperator(U*R_zm*S*U');

end
