#include<vector>
#include<stdlib.h>
#include<complex>
#include<iostream>
#include<cmath>

using namespace std;

extern "C" void zgeev_( char* jobvl, char* jobvr, int* n, complex<double>* a,
        int* lda, complex<double>* w, complex<double>* vl, int* ldvl, complex<double>* vr, int* ldvr,
        complex<double>* work, int* lwork, double* rwork, int* info );

extern "C" void zheev_(char *jobz, char *uplo, int *n, complex<double> *a,int *lda, 
        double *w, complex<double> *work, int *lwork, double *rwork, int *info);

extern "C" int zcopy_(int *n, complex<double> *zx, int *incx, 
        complex<double> *zy, int *incy);

extern "C" int zgemm_(char *transa, char *transb, int *m, int *n, int *k, 
        complex<double> *alpha, complex<double> *a, int *lda, 
        complex<double> *b, int *ldb, complex<double> *beta, 
        complex<double> *c, int *ldc);
typedef vector< vector<double> > matrixd;
typedef vector<double> vectord;
typedef vector< vector<complex<double> > > matrixc;
typedef vector<complex<double> > vectorc;
typedef complex<double> complexd;

#define PI 3.14159265
complexd I = complexd(0,1);

void CallZheev( matrixc& H, matrixc&, int m);

int  gcd(int  a,int  b){return b ? gcd(b, a%b) : a;}

void chop1( double&  v){if(abs(v)<1e-8)v=0;}

void chop1( complex<double>&  v){
    if(abs(v.real())<1e-8)v-=v.real();
    if(abs(v.imag())<1e-8)v-=I*v.imag();
}

void chop( vector<double > & v,int n){ 
    for(int i=0; i<n*n; i++)
        chop1(v[i]);
}

void chop( vector<complex<double> > & v,int n){ 
    for(int i=0; i<n*n; i++)
        chop1(v[i]);
}

void chop( double * v, int n){ 
    for(int i=0; i<n*n; i++)
        chop1(v[i]);
}

void chop( complex<double> * v, int n){ 
    for(int i=0; i<n*n; i++)
        chop1(v[i]);
}
void dbg(const char * c, complex<double> * v, int n){
    cout<<c;//chop(d);
    for (int i=0; i < n; i++) //Pass eigenvectors back, stored as ROWS !!!
        cout<<v[i]<<" ";
    cout<<endl;
}

inline double kron(int j, int k){return (double)(j==k);}

bool unique(int p, int q){
    return gcd(p, q)==1;
}

int main(){
    //cout<<"Complex"<<I<<endl;

    int q=751,p=1; int kx_max =1;
    double phi = (double) p / (double) q;
    for(q=2;q<60;q++){
        for(p=1;p<q;p++){
            if(unique(p,q)){
                for(int kx_inc=0;kx_inc<kx_max;kx_inc++){
                    double kx=0, ky=0;
                    double phi = (double) p / (double) q;
                    matrixc H(q,vectorc(q,0));
                    matrixc H2(q,vectorc(q,0));
                    for(int j=0;j<q; j++){
                        for(int i=0;i<q;i++){
                            H[i][j]		=kron(i,j)*2*cos(ky - 2*PI*(j+1)*phi);
                            H2[i][j]	=kron(i+1,j)+kron(i,j+1)+kron(i+q-1,j)*exp(-I*(q*kx))+kron(i,j+q-1)*exp(I*(q*kx));
                        }
                    }

                    CallZheev(H,H2,q);
                    for(int i=0;i<q;i++)
                        cout<<phi<<" "<<H[i][i].real()<<endl;//" "<<H[1][1]<<" "<<H[2][2]<<endl;
                }
            }
        }
    }
    //		cout<<H2[0][0]<<" "<<H2[0][1]<<" "<<H2[1][0]<<" "<<H2[1][1]<<endl;
    return 0;
}

void CallZheev( matrixc& H1, matrixc& H2, int m){    
    char jobz, uplo, conj='C', noop='N';
    int info, n=m, nn=n*n, incx=1, lwork = 3*n-1, lwork2 = 4*n, sizecomplex = sizeof(complex<double>);
    complex<double> alpha = 1, betaz = 0;
    jobz = 'V'; // V/N indicates that eigenvectors should/should not be calculated
    uplo = 'U'; // U/L indicated that the upper/lower triangle of the symmetric matrix is stored
    double *rwork = new double[3*n-2];
    double *rwork2 = new double[2*n];
    double *EN = new double[n];
    complex<double> *work = new complex<double>[lwork];
    complex<double> *work2 = new complex<double>[lwork2];
    complex<double> *a = new complex<double>[n*n];
    complex<double> *b = new complex<double>[n*n];
    complex<double> *c = new complex<double>[n*n];
    complex<double> *d = new complex<double>[n*n];
    complex<double> *e = new complex<double>[n*n];
    complex<double> *f = new complex<double>[n*n];
    complex<double> *dummy = new complex<double>[n*n];
    complex<double> *expw = new complex<double>[n];
    complex<double> * w = new complex<double>[n];	
    complex<double> * eigs = new complex<double>[n*n];	
    complex<double> *eomega = new complex<double>[n];

    for (int j=0; j < n; j++){ // transpose opration
        for (int i=0; i < n; i++){
            a[n*j+i] = H1[j][i];
            d[n*j+i] = H2[j][i];
        }
    }

    zheev_(&jobz, &uplo, &n, a, &n, EN, work, &lwork, rwork, &info);

    for (int i=0; i < n; i++) 
        for(int j = 0; j < n; j++)
            b[n*j+i] = exp(-I*EN[i]*PI*.5)* std::conj(a[n*i+j]);

    zgemm_(&noop,&noop, &n, &n, &n, &alpha, a, &n, b, &n, &betaz, c, &n);
    zheev_(&jobz, &uplo, &n, d, &n, EN, work, &lwork, rwork, &info);

    for (int i=0; i < n; i++) 
        for(int j = 0; j < n; j++)
            e[n*j+i] = exp(-I*EN[i]*PI*.5)*std::conj(d[n*i+j]);

    zgemm_(&noop,&noop,	&n, &n, &n, &alpha, d, &n, e, &n, &betaz, f, &n);
    zgemm_(&noop,&noop,	&n, &n, &n, &alpha, f, &n, c, &n, &betaz, a, &n);
    zgeev_( &noop, &jobz, &n, a, &n,  w, dummy, &n, eigs, &n,work2, &lwork2, rwork2, &info );

    for (int i=0; i < n; i++){ //Pass eigenvectors back, stored as ROWS !!!
        H1[i][i]=arg(w[i]);
        //				for(int j = 0; j < n; j++){
        //						H1[j][i] = c[n*i+j];
        //						H2[j][i] = f[n*i+j];
        //				}
    }
}

