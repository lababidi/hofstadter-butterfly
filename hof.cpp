#include<vector>
#include<stdlib.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<cmath>
#include<queue>

#ifdef __cplusplus
extern "C"
{
#endif
   #include <cblas.h>
#ifdef __cplusplus
}
#endif

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
#define N 3
typedef vector< vector<double> > matrixd;
typedef vector<double> vectord;
typedef vector< vector<complex<double> > > matrixc;
typedef vector<complex<double> > vectorc;
typedef complex<double> complexd;

#define PI 3.14159265
complexd I = complexd(0,1);

//void eigenvalues( matrixc& H, matrixc&, int m);
void eigenvalues( complex<double> a[],int n);

string calcH(int p, int q, int kn, ofstream &f);

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
    return (double) p/ (double) q<=.5 && gcd(p, q)==1 ;
}

int main(){

    ofstream f("but.csv");
    f<<"x,y"<<endl;
    #pragma omp parallel for
    for(int q=2;q<37;q++){
        string ff;
        for(int p=1;p<q;p++){
            if(unique(p,q)){
                ff+=calcH(p,q,20,f);
            }
        }
        f<<ff;
    }
    f.close();
    return 0;
}

struct xy{
    double phi;
    double energy;
    };

string calcH(int p, int q, int kn, ofstream &f){
    string ff = "";
    queue<xy> energies;
    ostringstream convert;
    double kr[kn];
    for(int kk=0; kk<kn; kk++)
        kr[kk]=acos(kk*2.0/kn-1)/q;

    for(int kyi=0; kyi<kn;kyi++){
        double ky = kr[kyi];
        for(int kxi=0; kxi<kyi; kxi++){
            double kx = kr[kxi];
            double phi = (double) p / (double) q;
            complex<double> *a = new complex<double>[q*q];
            for(int j=0;j<q; j++){
                for(int i=0;i<q;i++){
                    a[q*j+i]=kron(i,j)*2*cos(ky - 2*PI*(j+1)*phi)
                        +kron(i+1,j)+kron(i,j+1)
                        +kron(i+q-1,j)*exp(-I*(q*kx))
                        +kron(i,j+q-1)*exp(I*(q*kx));
                }
            }
            eigenvalues(a,q);
            for(int i=0;i<q;i++){
                energies.push((xy){phi,a[i*q+i].real()});
                convert<<phi<<","<<a[i*q+i].real()<<endl;
                convert<<(1-phi)<<","<<a[i*q+i].real()<<endl;
                ff+=convert.str();
            }
        }
    }
    return ff;
}




void eigenvalues( complex<double> a[],int n){    
    // V/N indicates that eigenvectors should/should not be calculated
    // U/L indicated that the upper/lower triangle of the symmetric matrix is stored
    char jobz='V', uplo='U', conj='C', noop='N';
    int info, lwork = 3*n-1;
    double *rwork = new double[3*n-2],
             *EN = new double[n];
    complex<double> *work = new complex<double>[lwork];

    zheev_(&jobz, &uplo, &n, a, &n, EN, work, &lwork, rwork, &info);
    for (int i=0; i < n; i++){ //Pass eigenvectors back, stored as ROWS !!!
        a[i*n+i]=EN[i];
    }
    return;
}


