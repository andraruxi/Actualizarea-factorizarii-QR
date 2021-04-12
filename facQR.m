function [ Q1,R1 ] = facQR( A, u, v )

    [Q,R] = qr(A);
    %[Qtest, Rtest] = qrupdate(Q, R, u, v);
    
    A = A + u * v';
    w = Q' * u;
    
    [m, n] = size(A);%am aflat dimensiunile matricei
    
    for k = m-1:-1:1%aduce vectorul w astfel il face sa aiba un sigur element pe primul rand
        J = eye(m);
        r = sign(w(k)) * sqrt(w(k)*w(k)+w(k+1)*w(k+1));
        c = w(k)/r;
        s = w(k+1)/r;%am folosit un algoritm putin mai complex pt rotatii deoarece vrem sa rotim toata matricea
        J(k,k) = c;%si nu doar mijlocul s,c
        J(k,k+1) = -s;
        J(k+1,k) = s;
        J(k+1,k+1) = c;
        
        w = J' * w;%foloseste matricele gives ca sa roteasca w 
                   %vrem sa il facem pe R superior Hessenberg de aceea
                   %folosim  rotatiile givens
        R = J' * R;
        Q = Q * J;
        %il inmultim acum pe Q ca sa tinem pasul si sa nu uitam sa il
        %inmultim la final
    end
    H1 = R + w * v';

    for k = 1:n%il facem pe H1 superior triunghiularea=>face givens doar pentru subdiagonala diagonalei principale(matrice Hessenberg) deoarece in rest avem zero=>mai eficient
        r = sign(H1(k,k))*sqrt(H1(k,k)*H1(k,k)+H1(k+1,k)*H1(k+1,k));
        c = H1(k,k) / r;
        s = H1(k+1,k) / r;
            
        G = eye(m);%rotatii givens pentru elementele de pe subdiagonala 
        G(k,k) = c;%la fel folosim un algoritm putin mai complex decat cel din curs
        G(k,k+1) = -s;
        G(k+1,k) = s;%facem rotatii k si k+1 cum ne sugereaza in culegere
        G(k+1,k+1) = c;
        
        H1 = G' * H1;
        Q = Q * G;
        
    end
    
    R1 = H1;%pe Q il inmultim cu toti givensii posibili
    Q1 = Q;
    
   
    
    
end
