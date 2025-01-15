
B= 200000000.0
B = 2.000015228986740e+08

f_c = 300000000.0
f_c = 3.000007614027709e+08
p=2

2*pi*f_c./(B.*p)
% outside debuggin we get
4.712388980384690
%inside debuggin it gives
4.712365058198953


% its beciase of slightly different values of B and fc


3.000007614027709e+08

% in python 
%2*np.pi*f_c/(B*p)
% gives
%4.71238898038469
