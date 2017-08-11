#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;
int main(int argc, char *argv[]) {
	
	int ratio = 6;
	double small_tau = 0.2/double(ratio);
	int Step = 886006;
	
	double cu_time = 177201.23333;
	
	cout<<"Before taking the floor function: "<<((cu_time - small_tau * ratio * Step)/small_tau)<<endl;
	cout<<"After taking floor function: "<<int((cu_time - small_tau * ratio * Step)/small_tau)<<endl;
	
	bool found = false;
	for(int i = 0 ; i < 10 ; i ++) {
		for(int j = 0 ; j < 10 ; j++) {
			if(i==0 && j==0) {
				found = true;
				break;
			}
		}
		if(found) break;
		cout<<"OKOKOKO";
	}
	//It should be in bucket 1 in my opinion.
}