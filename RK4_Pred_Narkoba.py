import matplotlib.pyplot as plt

#  KESTABILAN DINAMIK KELOMPOK PECANDU NARKOBA DENGAN MELIBATKAN TREATMENT PENCEGAHAN EDUKASI, REHABILITASI DAN PENGAWASAN INTERAKSI 
#  MENGGUNAKAN METODE RK4

# Fungsi SPD Tanpa Control
def FN(S,P,a,b,c,d,f):
    return (a-(b*S)-(c*S)-((f*S*P)/367602)+(d*P))

def GN(S,P,Pa,b,d,f,z,i,j):
    return ((f*S*P)/367602)+((z*P*Pa)/367602)-((b+j)*P)-(i*P)-(d*P)

def HN(S,P,Sa,b,c,g):
    return (c*S)-(b*Sa)-((g*Sa*P)/367602)

def IN(P,Sa,Pa,b,e,g,z,i):
    return i*P-((z*P*Pa)/367602)-b*Pa-e*Pa+((g*Sa*P)/367602)

def JN(Pa,Ra,b,e):
    return (Pa*e)-(b*Ra)

# Fungsi State SPD Control Terhadap C1, C2 dan K
def FB(S,P,a,b,c,d,f,C1):
    return (a-(b*S)-(c*C1*S)-((f*S*P)/367602)+(d*P))

def GB(S,P,Pa,b,d,f,z,i,j,C2,K):
    return ((f*S*P)/367602)+((z*K*P*Pa)/367602)-((b+j)*P)-(i*C2*P)-(d*P)

def HB(S,P,Sa,b,c,g,C1):
    return (c*C1*S)-(b*Sa)-((g*Sa*P)/367602)

def IB(P,Sa,Pa,b,e,g,z,i,C2,K):
    return (i*C2*P)-((z*K*P*Pa)/367602)-(b*Pa)-(e*Pa)+((g*Sa*P)/367602)

def JB(Pa,Ra,b,e):
    return (Pa*e)-(b*Ra)

# Fungsi Co-State SPD Control Terhadap C1, C2 dan K
def FB2(lambda1,lambda2,lambda3,P,b,c,f,C1):
    return (((c*C1+b)*lambda1-lambda3*c*C1)+((P*f*(lambda1-lambda2))/367602))

def GB2(lambda1,lambda2,lambda3,lambda4,S,Sa,Pa,b,d,f,g,z,i,j,C2,K):
    return ((((i*C2+b+d+j)*lambda2-i*C2*lambda4-d*lambda1-1)*367602)+(-K*Pa*z-S*f)*lambda2
            +(K*Pa*z-Sa*g)*lambda4+S*lambda1*f+lambda3*g*Sa)/367602

def HB2(lambda3,lambda4,P,b,g):
    return (lambda3*(367602*b+P*g)/367602)-(lambda4*g*P/367602)

def IB2(lambda2,lambda4,lambda5,P,b,e,z,K):
    return -((lambda2*z*K*P)/367602)-lambda4*((-z*K*P/367602)-b-e)-lambda5*e

def JB2(lambda5,b):
    return lambda5*b
 
# NILAI PERAMETER PADA SPD
a0 = 3.857
b0 = 2.356
c0 = 0.142
d0 = 1-0.672
e0 = 0.672
f0 = 0.572
g0 = 0.427
z0 = 0.011
i0 = 0.712
j0 = 0.031

# DATA AWAL SUBPOPULASI PADA SPD KONTROL TERHADAP C1, C2 dan K
S0_B = 293943
P0_B = 300
Pa0_B = 1550
Sa0_B = 72127
Ra0_B = 1530

# DATA AWAL SUBPOPULASI PADA SPD TANPA PENGONTROL
S0_N = 293943
P0_N = 300
Pa0_N = 1550
Sa0_N = 72127
Ra0_N = 1530

# NILAI AWAL YANG DIBERIKAN TERHADAP PENGONTROL
C10_B = 0.1
C20_B = 0.1
K0_B = 0.1

# BATAS NILAI MINIMUM DAN MAXIMUM PADA PENGONTROL
C1max = 0.9
C1min = 0.1
C2max = 0.9
C2min = 0.1
Kmax = 0.9
Kmin = 0.1

# NILAI AWAL LAMBDA
lambda10_B = 0
lambda20_B = 0
lambda30_B = 0
lambda40_B = 0
lambda50_B = 0

# PREDIKSI PECANDU NARKOBA PADA HARI Ke-n
T = 1
n = 365*10   
h = 1/(365*5)
h2 = h/2

# UPDATE DATA SUBPOPULASI SPD TANPA KONTROL TERHADAP WAKTU-t / dimana kondisi ini t = n
data_S_N = []
data_P_N = []
data_Sa_N = []
data_Pa_N = []
data_Ra_N = []

# UPDATE DATA SUBPOPULASI SPD KONTROL TERHADAP C1, C2 dan K TERHADAP WAKTU-t / dimana konidis ini t = n
data_S_B = []
data_P_B = []
data_Sa_B = []
data_Pa_B = []
data_Ra_B = []

# UPDATE DATA WAKTU YANG BERJALAN DARI 0 < t < n
data_waktu = []

# WAKTU AWAL DIMULAI
i = 0


def rk4(C10_B,C20_B,K0_B,C1max,C1min,C2max,C2min,lambda10_B,lambda20_B,lambda30_B,lambda40_B,lambda50_B,S0_B,P0_B,Sa0_B,Pa0_B,Ra0_B,S0_N,P0_N,Sa0_N,Pa0_N,Ra0_N,n,
        data_S_N,data_P_N,data_Sa_N,data_Pa_N,data_Ra_N, data_S_B,data_P_B,data_Sa_B,data_Pa_B,data_Ra_B, data_waktu,t0):


    print('\n---------------------------------------------Solution Of Variabel No-Control----------------------------------------------')
    print('---------------------------------------------------------------------------------------------------------------------------') 
    print('t(0)\tS(0)\t\tP(0)\tSa(0)\t\tPa(0)\tRa(0)\tS(n)\t\tP(n)\tSa(n)\t\tPa(n)\tRa(n)')
    print('---------------------------------------------------------------------------------------------------------------------------')

    for i in range(n):
            data_S_N += [S0_N]
            data_P_N += [P0_N]
            data_Sa_N += [Sa0_N]
            data_Pa_N += [Pa0_N]
            data_Ra_N += [Ra0_N]
            data_waktu += [t0]

            k1x1 =  (FN(S0_N,P0_N,a0,b0,c0,d0,f0))
            k1x2 =  (GN(S0_N,P0_N,Pa0_N,b0,d0,f0,z0,i0,j0))
            k1x3 =  (HN(S0_N,P0_N,Sa0_N,b0,c0,g0))
            k1x4 =  (IN(P0_N,Sa0_N,Pa0_N,b0,e0,g0,z0,i0))
            k1x5 =  (JN(Pa0_N,Ra0_N,b0,e0))

            k2x1 =  (FN(S0_N+((h2/2)*k1x1),P0_N+((h2/2)*k1x1),a0/2,b0/2,c0/2,d0/2,f0/2))
            k2x2 =  (GN(S0_N+((h2/2)*k1x2),P0_N+((h2/2)*k1x2),Pa0_N+((h2/2)*k1x2),b0/2,d0/2,f0/2,z0/2,i0/2,j0/2))
            k2x3 =  (HN(S0_N+((h2/2)*k1x3),P0_N+((h2/2)*k1x3),Sa0_N+((h2/2)*k1x3),b0/2,c0/2,g0/2))
            k2x4 =  (IN(P0_N+((h2/2)*k1x4),Sa0_N+((h2/2)*k1x4),Pa0_N+((h2/2)*k1x4),b0/2,e0/2,g0/2,z0/2,i0/2))
            k2x5 =  (JN(Pa0_N+((h2/2)*k1x5),Ra0_N+((h2/2)*k1x5),b0/2,e0/2))

            k3x1 =  (FN(S0_N+((h2/2)*k2x1),P0_N+((h2/2)*k2x1),a0/2,b0/2,c0/2,d0/2,f0/2))
            k3x2 =  (GN(S0_N+((h2/2)*k2x2),P0_N+((h2/2)*k2x2),Pa0_N+((h2/2)*k2x2),b0/2,d0/2,f0/2,z0/2,i0/2,j0/2))
            k3x3 =  (HN(S0_N+((h2/2)*k2x3),P0_N+((h2/2)*k2x3),Sa0_N+((h2/2)*k2x3),b0/2,c0/2,g0/2))
            k3x4 =  (IN(P0_N+((h2/2)*k2x4),Sa0_N+((h2/2)*k2x4),Pa0_N+((h2/2)*k2x4),b0/2,e0/2,g0/2,z0/2,i0/2))
            k3x5 =  (JN(Pa0_N+((h2/2)*k2x5),Ra0_N+((h2/2)*k2x5),b0/2,e0/2))

            k4x1 =  (FN(S0_N+(h2*k3x1),P0_N+(h2*k3x1),a0,b0,c0,d0,f0))
            k4x2 =  (GN(S0_N+(h2*k3x2),P0_N+(h2*k3x2),Pa0_N+(h2*k3x2),b0,d0,f0,z0,i0,j0))
            k4x3 =  (HN(S0_N+(h2*k3x3),P0_N+(h2*k3x3),Sa0_N+(h2*k3x3),b0,c0,g0))
            k4x4 =  (IN(P0_N+(h2*k3x4),Sa0_N+(h2*k3x4),Pa0_N+(h2*k3x4),b0,e0,g0,z0,i0))
            k4x5 =  (JN(Pa0_N+(h*k3x5),Ra0_N+(h*k3x5),b0,e0))

            k1 = h*(k1x1+2*k2x1+2*k3x1+k4x1)/6
            k2 = h*(k1x2+2*k2x2+2*k3x2+k4x2)/6
            k3 = h*(k1x3+2*k2x3+2*k3x3+k4x3)/6
            k4 = h*(k1x4+2*k2x4+2*k3x4+k4x4)/6
            k5 = h*(k1x5+2*k2x5+2*k3x5+k4x5)/6
            
            Sn_N = S0_N + k1
            Pn_N = P0_N + k2
            San_N = Sa0_N + k3
            Pan_N = Pa0_N + k4
            Ran_N = Ra0_N + k5
            
            P0_N = Pn_N
            S0_N = Sn_N
            Sa0_N = San_N
            Pa0_N = Pan_N
            Ra0_N = Ran_N

            t0 = i

            print('%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'%(t0,S0_N,P0_N,Sa0_N,Pa0_N,Ra0_N,Sn_N,Pn_N,San_N,Pan_N,Ran_N))
            print('---------------------------------------------------------------------------------------------------------------------------')
    
    print('\nAt S=%.2f, P=%.2f, Sa=%.2f, Pa=%.2f, Ra=%.2f'%(Sn_N,Pn_N,San_N,Pan_N,Ran_N))
    print('\nAt Waktu = %i'%t0)
    
    print('\n----------------------------------------------Solution Of Variabel C1, C2 and K-------------------------------------------')
    print('---------------------------------------------------------------------------------------------------------------------------') 
    print('t(0)\tS(0)\t\tP(0)\tSa(0)\t\tPa(0)\tRa(0)\tS(n)\t\tP(n)\tSa(n)\t\tPa(n)\tRa(n)\tC1(n)\tC2(n)\tK(t)')
    print('---------------------------------------------------------------------------------------------------------------------------')

    for i in range(n):
            data_S_B += [S0_B]
            data_P_B += [P0_B]
            data_Sa_B += [Sa0_B]
            data_Pa_B += [Pa0_B]
            data_Ra_B += [Ra0_B]

            k1x1 =  (FB(S0_B,P0_B,a0,b0,c0,d0,f0,C10_B))
            k1x2 =  (GB(S0_B,P0_B,Pa0_B,b0,d0,f0,z0,i0,j0,C20_B,K0_B))
            k1x3 =  (HB(S0_B,P0_B,Sa0_B,b0,c0,g0,C10_B))
            k1x4 =  (IB(P0_B,Sa0_B,Pa0_B,b0,e0,g0,z0,i0,C20_B,K0_B))
            k1x5 =  (JB(Pa0_B,Ra0_B,b0,e0))

            k2x1 =  (FB(S0_B+((h2/2)*k1x1),P0_B+((h2/2)*k1x1),a0/2,b0/2,c0/2,d0/2,f0/2,C10_B/2))
            k2x2 =  (GB(S0_B+((h2/2)*k1x2),P0_B+((h2/2)*k1x2),Pa0_B+((h2/2)*k1x2),b0/2,d0/2,f0/2,z0/2,i0/2,j0/2,C20_B/2,K0_B/2))
            k2x3 =  (HB(S0_B+((h2/2)*k1x3),P0_B+((h2/2)*k1x3),Sa0_B+((h2/2)*k1x3),b0/2,c0/2,g0/2,C10_B/2))
            k2x4 =  (IB(P0_B+((h2/2)*k1x4),Sa0_B+((h2/2)*k1x4),Pa0_B+((h2/2)*k1x4),b0/2,e0/2,g0/2,z0/2,i0/2,C20_B/2,K0_B/2))
            k2x5 =  (JB(Pa0_B+((h2/2)*k1x5),Ra0_B+((h2/2)*k1x5),b0/2,e0/2))

            k3x1 =  (FB(S0_B+((h2/2)*k2x1),P0_B+((h2/2)*k2x1),a0/2,b0/2,c0/2,d0/2,f0/2,C10_B/2))
            k3x2 =  (GB(S0_B+((h2/2)*k2x2),P0_B+((h2/2)*k2x2),Pa0_B+((h2/2)*k2x2),b0/2,d0/2,f0/2,z0/2,i0/2,j0/2,C20_B/2,K0_B/2))
            k3x3 =  (HB(S0_B+((h2/2)*k2x3),P0_B+((h2/2)*k2x3),Sa0_B+((h2/2)*k2x3),b0/2,c0/2,g0/2,C10_B/2))
            k3x4 =  (IB(P0_B+((h2/2)*k2x4),Sa0_B+((h2/2)*k2x4),Pa0_B+((h2/2)*k2x4),b0/2,e0/2,g0/2,z0/2,i0/2,C20_B/2,K0_B/2))
            k3x5 =  (JB(Pa0_B+((h2/2)*k2x5),Ra0_B+((h2/2)*k2x5),b0/2,e0/2))

            k4x1 =  (FB(S0_B+(h2*k3x1),P0_B+(h2*k3x1),a0,b0,c0,d0,f0,C10_B))
            k4x2 =  (GB(S0_B+(h2*k3x2),P0_B+(h2*k3x2),Pa0_B+(h2*k3x2),b0,d0,f0,z0,i0,j0,C20_B,K0_B))
            k4x3 =  (HB(S0_B+(h2*k3x3),P0_B+(h2*k3x3),Sa0_B+(h2*k3x3),b0,c0,g0,C10_B))
            k4x4 =  (IB(P0_B+(h2*k3x4),Sa0_B+(h2*k3x4),Pa0_B+(h2*k3x4),b0,e0,g0,z0,i0,C20_B,K0_B))
            k4x5 =  (JB(Pa0_B+(h*k3x5),Ra0_B+(h*k3x5),b0,e0))

            k1 = h*(k1x1+2*k2x1+2*k3x1+k4x1)/6
            k2 = h*(k1x2+2*k2x2+2*k3x2+k4x2)/6
            k3 = h*(k1x3+2*k2x3+2*k3x3+k4x3)/6
            k4 = h*(k1x4+2*k2x4+2*k3x4+k4x4)/6
            k5 = h*(k1x5+2*k2x5+2*k3x5+k4x5)/6
            
            Sn_B = S0_B + k1
            Pn_B = P0_B + k2
            San_B = Sa0_B + k3
            Pan_B = Pa0_B + k4
            Ran_B = Ra0_B + k5
            
            P0_B = Pn_B
            S0_B = Sn_B
            Sa0_B = San_B
            Pa0_B = Pan_B
            Ra0_B = Ran_B

            w1y1 =  FB2(lambda10_B,lambda20_B,lambda30_B,P0_B,b0,c0,f0,C10_B)
            w1y2 =  GB2(lambda10_B,lambda20_B,lambda30_B,lambda40_B,S0_B,Sa0_B,Pa0_B,b0,d0,f0,g0,z0,i0,j0,C20_B,K0_B)
            w1y3 =  HB2(lambda30_B,lambda40_B,P0_B,b0,g0)
            w1y4 =  IB2(lambda20_B,lambda40_B,lambda50_B,P0_B,b0,e0,z0,K0_B)
            w1y5 =  JB2(lambda50_B,b0)

            w2y1 =  FB2(lambda10_B-h2*w1y1,lambda20_B-h2*w1y1,lambda30_B-h2*w1y1,P0_B/2,b0/2,c0/2,f0/2,C10_B/2)
            w2y2 =  GB2(lambda10_B-h2*w1y2,lambda20_B-h2*w1y2,lambda30_B-h2*w1y2,lambda40_B-h2*w1y2,S0_B/2,Sa0_B/2,Pa0_B/2,b0/2,d0/2,f0/2,g0/2,z0/2,i0/2,j0/2,C20_B/2,K0_B/2)
            w2y3 =  HB2(lambda30_B-h2*w1y3,lambda40_B-h2*w1y3,P0_B/2,b0/2,g0/2)
            w2y4 =  IB2(lambda20_B-h2*w1y4,lambda40_B-h2*w1y4,lambda50_B-h2*w1y4,P0_B/2,b0/2,e0/2,z0/2,K0_B/2)
            w2y5 =  JB2(lambda50_B-h2*w1y5,b0/2)

            w3y1 =  FB2(lambda10_B-h2*w2y1,lambda20_B-h2*w2y1,lambda30_B-h2*w2y1,P0_B/2,b0/2,c0/2,f0/2,C10_B/2)
            w3y2 =  GB2(lambda10_B-h2*w2y2,lambda20_B-h2*w2y2,lambda30_B-h2*w2y2,lambda40_B-h2*w2y2,S0_B/2,Sa0_B/2,Pa0_B/2,b0/2,d0/2,f0/2,g0/2,z0/2,i0/2,j0/2,C20_B/2,K0_B/2)
            w3y3 =  HB2(lambda30_B-h2*w2y3,lambda40_B-h2*w2y3,P0_B/2,b0/2,g0/2)
            w3y4 =  IB2(lambda20_B-h2*w2y4,lambda40_B-h2*w2y4,lambda50_B-h2*w2y4,P0_B/2,b0/2,e0/2,z0/2,K0_B/2)
            w3y5 =  JB2(lambda50_B-h2*w2y5,b0/2)

            w4y1 =  FB2(lambda10_B-w3y1,lambda20_B-w3y1,lambda30_B-w3y1,P0_B,b0,c0,f0,C10_B)
            w4y2 =  GB2(lambda10_B-w3y2,lambda20_B-w3y2,lambda30_B-w3y2,lambda40_B-w3y2,S0_B,Sa0_B,Pa0_B,b0,d0,f0,g0,z0,i0,j0,C20_B,K0_B)
            w4y3 =  HB2(lambda30_B-w3y3,lambda40_B-w3y3,P0_B,b0,g0)
            w4y4 =  IB2(lambda20_B-w3y4,lambda40_B-w3y4,lambda50_B-w3y4,P0_B,b0,e0,z0,K0_B)
            w4y5 =  JB2(lambda50_B-w3y5,b0)

            w1 = h*(w1y1+2*w2y1+2*w3y1+w4y1)/6
            w2 = h*(w1y2+2*w2y2+2*w3y2+w4y2)/6
            w3 = h*(w1y3+2*w2y3+2*w3y3+w4y3)/6
            w4 = h*(w1y4+2*w2y4+2*w3y4+w4y4)/6
            w5 = h*(w1y5+2*w2y5+2*w3y5+w4y5)/6

            lambda1n_B = lambda10_B - w1
            lambda2n_B = lambda20_B - w2
            lambda3n_B = lambda30_B - w3
            lambda4n_B = lambda40_B - w4
            lambda5n_B = lambda50_B - w5

            lambda1n_B = max(C1min,min(lambda10_B,C1max))
            lambda2n_B = max(C1min,min(lambda20_B,C1max))
            lambda3n_B = max(C1min,min(lambda30_B,C1max))
            lambda4n_B = max(C1min,min(lambda40_B,C1max))
            lambda5n_B = max(C1min,min(lambda50_B,C1max))

            C1n_B = min(C1max,max(C1min,(41739.906*lambda1n_B)-(41739.906*lambda3n_B)))
            C2n_B = min(C2max,max(C2min,(71.200*lambda2n_B)-(71.200*lambda5n_B)))
            Kn_B = min(0.9,max(0.1,(-0.4638168453e-2*lambda2n_B)+(0.4638168453e-2*lambda5n_B)))

            C10_B = C10_B + C1n_B
            C20_B = C20_B + C2n_B
            K0_B = K0_B + Kn_B

            lambda10_B = lambda1n_B
            lambda20_B = lambda2n_B
            lambda30_B = lambda3n_B
            lambda40_B = lambda4n_B
            lambda50_B = lambda5n_B

            t0 = i
            

            print('%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.7f\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'%(t0,S0_B,P0_B,Sa0_B,Pa0_B,Ra0_B,Sn_B,Pn_B,San_B,Pan_B,Ran_B,C1n_B,C2n_B,Kn_B))
            print('---------------------------------------------------------------------------------------------------------------------------')


    print('\nAt S=%.2f, P=%.2f, Sa=%.2f, Pa=%.2f, Ra=%.2f'%(Sn_B,Pn_B,San_B,Pan_B,Ran_B))
    print('\nAt L1=%.2f, L2=%.2f, L3=%.2f, L4=%.2f, L5=%.2f' %(lambda10_B,lambda20_B,lambda30_B,lambda40_B,lambda50_B))
    print('\nAt C1(0)=%.2f, C2(0)=%.2f'%(C10_B, C20_B))
    print('\nAt K = %.2f'%(K0_B))
    print('\nAt Waktu = %i'%t0)

rk4(C10_B,C20_B,K0_B,C1max,C1min,C2max,C2min,lambda10_B,lambda20_B,lambda30_B,lambda40_B,lambda50_B,S0_B,P0_B,Sa0_B,Pa0_B,Ra0_B,S0_N,P0_N,Sa0_N,Pa0_N,Ra0_N,n,
    data_S_N,data_P_N,data_Sa_N,data_Pa_N,data_Ra_N,data_S_B,data_P_B,data_Sa_B,data_Pa_B,data_Ra_B, data_waktu, i)


# Ploting Grafik SPD Tanpa Kontrol dan Grafik SPD Dengan Pengontrol C1, C2 dan K "SUBPOPULASI PECANDU NARKOBA"
plt.plot(data_waktu,data_P_N, color='red', linestyle='dotted',  linewidth = 2,label = 'Tanpa Pengontrol')
plt.plot(data_waktu,data_P_B, color='black', linestyle='dotted',  linewidth = 2,label ='Dengan Pengontrol C1, C2 dan K')

# Batas Garis x
plt.xlim(0,2000)

# Batas Garis y
plt.ylim(0,400)

# Penamaan Pada Garis x
plt.xlabel('Waktu (Hari)')

# Penamaan Pada Garis y
plt.ylabel('Populasi Masyarakat')

# Judul Pada Gambar Grafik
plt.title('Subpopulasi Pecandu Narkoba (P)')

# Keterangan Garis 
plt.legend()

# Menampilkan Hasil Ploting
plt.show()