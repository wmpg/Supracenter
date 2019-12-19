#include <iostream>
#include <iterator>
#include <math.h>

#define PI 3.14159265359
#define PI_2 1.57079632679

int cscan(float S[3], float D[3], float z_profile[][4], \
	int z_len, bool wind=true, int n_theta=90, int n_phi=90, float h_tol=1e-5, float v_tol=1000){

	if (!wind) {
		// std::cout<< "no wind";
		for (int i=0; i < z_len; i++) {

			z_profile[i][2] = 0;
		}
	}
	
	float d_theta = PI/n_theta;
	float d_phi   = PI_2/n_phi;

	float dx = D[0] - S[0];
	float dy = D[1] - S[1];



	return 0;
}



int main(){
	float S[3] = {0, 0, 5.5};
	float D[3] = {0, 0, 5.5};
	int z_len = 2;
	float z_profile[z_len][4] = {{0, 0, 0, 0}, {1, 1, 1.5, 1}};
	cscan(S, D, z_profile, z_len, false);
}