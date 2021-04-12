m=8;
n=7;
A=randn(m,n);
u=randn(m,1);
v=randn(n,1);
[Q,R]=qr(A);
[Qt,Rt]=qrupdate(Q,R,u,v);
[Q1,R1]=facQR(A,u,v);
[Q,R]=qr(A+u*v');
%Q1 da identic cu Qt 
%cand verificam Q cu Q1  avem uneori insa semne diferite!
Q
Q1


