v_c = input('clay content=');
I = eye (6,6);
alpha1 = 1;
alpha2 = 1/1000;
k_q1 = 37;
miu_q1 = 44;
rho_q = 2.65;
k_c1 = 22;
miu_c1 = 10;
rho_c = 2.5;
M = 100;
N = 100;
for m = 1:M
    phi(m) = 0.2*m/M;
    for n = 1:N
        eps(n) = 0.1*(n/N);
        v2 = (4/3)*pi*eps(n)*alpha2;
        v_q = 1-v_c;
        v_1 = v_q/(1-phi(m)-v2);
        v_2 = v_c/(1-phi(m)-v2);
        M_a = (k_q1*k_c1)/((k_c1*v_1)+(k_q1*v_2));
        M_b = (v_1*k_q1)+(v_2*k_c1);
        k_m = (M_a + M_b)/2;
        M_a1 = (miu_q1*miu_c1)/((miu_c1*v_1)+(miu_q1*v_2));
        M_b1 = (v_1*miu_q1)+(v_2*miu_c1);
        miu_m = (M_a1 + M_b1)/2;
        Lam = k_m - (2/3*miu_m);
        c11= Lam + 2*miu_m;
        c12 = Lam;
        c44 = miu_m;
        C_0 = [
    c11 c12 c12 0 0 0;
    c12 c11 c12 0 0 0;
    c12 c12 c11 0 0 0;
    0 0 0 2*c44 0 0;
    0 0 0 0 2*c44 0;
    0 0 0 0 0 2*c44
    ];
S_0 = inv (C_0);
G1 = Gtensor(C_0, 1/alpha1);
G2 = Gtensor(C_0, 1/alpha2);
k1 = inv(I + G1*C_0)*S_0;
k2 = inv(I + G2*C_0)*S_0;
k2 = iso_av(k2);
s_eff = S_0 + phi(m)*k1 + v2*k2;
c_eff = inv (s_eff);
rho_eff = (rho_q*(1-phi(m) - v_c) + v_c*rho_c)/(1 - phi(m));
c_eff(4,4) = (c_eff(4,4))/2;
Vp (m,n) = sqrt((c_eff(1,1))/(rho_eff));
Vs (m,n) = sqrt((c_eff(4,4)/2)/(rho_eff));
    end
end
figure(1);
surf(phi,eps,Vp);
xlabel('porosity');
ylabel('Crack density');
zlabel('P-wave velocity (Km/sec)');
figure(2);
surf(phi,eps,Vs);
xlabel('porosity');
ylabel('Crack density');
zlabel('S-wave velocity (Km/sec)');