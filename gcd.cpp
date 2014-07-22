#include<iostream>

using namespace std;

int  gcd(int  a,int  b, int & steps){return b ? gcd(b, a%b,++steps) : a;}

int main(){
	int steps = 1;
	for(int x=1; x<400; x++){
			for (int y=1; y<400; y++){
					steps = 1;
					cout<<gcd(y,x,steps)<<" "<<x<<" "<<y<<" "<<steps<<endl;
			}	
			cout<<endl;
	}

	return 0;
}
