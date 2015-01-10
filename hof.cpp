#include<vector>
#include<stdlib.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<cmath>
#include<queue>

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

//void CallZheev( matrixc& H, matrixc&, int m);
void CallZheev( complex<double> a[],int n);

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
                //if(p==1 || p==q-1 || p==(q-1)/2) calcH(p,q,16,f);
                //else calcH(p,q,1,f);
                ff+=calcH(p,q,20,f);
            }
        }
        f<<ff;
    }
    f.close();
    //		cout<<H2[0][0]<<" "<<H2[0][1]<<" "<<H2[1][0]<<" "<<H2[1][1]<<endl;
    return 0;
}

string calcH(int p, int q, int kn, ofstream &f){
    string ff = "";
    ostringstream convert;
    double kr[kn];
    for(int kk=0; kk<kn; kk++)
        kr[kk]=acos(kk*2.0/kn-1)/q;

    for(int kyi=0; kyi<kn;kyi++){
        double ky = kr[kyi];
        for(int kxi=0; kxi<kyi; kxi++){
            double kx = kr[kxi];
            double phi = (double) p / (double) q;
            matrixc H(q,vectorc(q,0));
            matrixc H2(q,vectorc(q,0));
            complex<double> *a = new complex<double>[q*q];
            for(int j=0;j<q; j++){
                for(int i=0;i<q;i++){
                    /*H[i][j]*/a[q*j+i]		=kron(i,j)*2*cos(ky - 2*PI*(j+1)*phi)
                        +kron(i+1,j)+kron(i,j+1)+kron(i+q-1,j)*exp(-I*(q*kx))+kron(i,j+q-1)*exp(I*(q*kx));
                }
            }
            CallZheev(a,q);
            for(int i=0;i<q;i++){
                convert<<phi<<","<<a[i*q+i].real()<<endl;
                convert<<(1-phi)<<","<<a[i*q+i].real()<<endl;
                ff+=convert.str();
            }
        }
    }
    return ff;
}




//void CallZheev( matrixc& H1, matrixc& H2, int n){    
void CallZheev( complex<double> a[],int n){    
    char jobz, uplo, conj='C', noop='N';
    int info, nn=n*n, incx=1, lwork = 3*n-1, lwork2 = 4*n, sizecomplex = sizeof(complex<double>);
    complex<double> alpha = 1, betaz = 0;
    jobz = 'V'; // V/N indicates that eigenvectors should/should not be calculated
    uplo = 'U'; // U/L indicated that the upper/lower triangle of the symmetric matrix is stored
    double *rwork = new double[3*n-2];
    double *rwork2 = new double[2*n];
    double *EN = new double[n];
    complex<double> *work = new complex<double>[lwork];
    complex<double> *work2 = new complex<double>[lwork2];
    //complex<double> *a = new complex<double>[n*n];
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

    /*
       for (int j=0; j < n; j++){ // transpose opration
       for (int i=0; i < n; i++){
       a[n*j+i] = H1[j][i];
       d[n*j+i] = H2[j][i];
       }
       }*/

    zheev_(&jobz, &uplo, &n, a, &n, EN, work, &lwork, rwork, &info);
    for (int i=0; i < n; i++){ //Pass eigenvectors back, stored as ROWS !!!
        //H1[i][i]=EN[i];
        a[i*n+i]=EN[i];
    }
    return;
}


