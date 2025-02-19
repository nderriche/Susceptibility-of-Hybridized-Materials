#include <stdio.h>
#include <math.h>


////////
// 1s function norm
double norm_1s(int n, double *x, void *imported_data){
	double a_0 = *((double *)(imported_data)+3);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	return pow(1/(sqrt(M_PI)*pow(a_0,1.5)),2) * exp(-2*r/a_0);
}


// 2s function norm
double norm_2s(int n, double *x, void *imported_data){
	double a_0 = *((double *)(imported_data)+3);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	return pow((2-r/a_0)/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * exp(-r/a_0);
}
////////


// 1s-1s
double real_1s1s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(sqrt(M_PI)*pow(a_0,1.5)),2) * exp(-r/(a_0)) * exp(-r_neighbor/(a_0)) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double imag_1s1s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(sqrt(M_PI)*pow(a_0,1.5)),2) * exp(-r/(a_0)) * exp(-r_neighbor/(a_0)) * sin(qx*x[0] + qy*x[1] + qz*x[2]);
}



// 2s-2s
double real_2s2s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r/a_0) * (2-r_neighbor/a_0) * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double imag_2s2s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r/a_0) * (2-r_neighbor/a_0) * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * sin(qx*x[0] + qy*x[1] + qz*x[2]);
}




// 2p-2p (Aligned) 
double real_2p2p_aligned(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (x[2]/a_0) * ((x[2]-az)/a_0) * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double imag_2p2p_aligned(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (x[2]/a_0) * ((x[2]-az)/a_0) * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * sin(qx*x[0] + qy*x[1] + qz*x[2]);

}



// 2s-2p 

double real_2s2p(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r_neighbor/a_0) * x[2]/a_0 * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double imag_2s2p(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r_neighbor/a_0) * x[2]/a_0 * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * sin(qx*x[0] + qy*x[1] + qz*x[2]);
}



// 2p-2s 

double real_2p2s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r/a_0) * (x[2]-az)/a_0 * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double imag_2p2s(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);

	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double r_neighbor = sqrt((x[0] - ax)*(x[0] - ax) + (x[1] - ay)*(x[1] - ay) + (x[2] - az)*(x[2] - az));

	return pow(1/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * (2-r/a_0) * (x[2]-az)/a_0 * exp(-r/(2*a_0)) * exp(-r_neighbor/(2*a_0)) * sin(qx*x[0] + qy*x[1] + qz*x[2]);
}




/*
// 2p-2p (Orthogonal) Onsite
double onsite_2p2p_orthogonal_real(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	return pow((x[2]/a_0)/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * exp(-r/a_0) * cos(qx*x[0] + qy*x[1] + qz*x[2]);
}

double onsite_2p2p_orthogonal_imag(int n, double *x, void *imported_data){
	double qx =  *(double *)imported_data, qy = *((double *)(imported_data)+1), qz = *((double *)(imported_data)+2), a_0 = *((double *)(imported_data)+3);
	double ax = *((double *)(imported_data)+4), ay = *((double *)(imported_data)+5), az = *((double *)(imported_data)+6);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	return pow((x[2]/a_0)/(4.0*sqrt(2*M_PI)*pow(a_0,1.5)),2) * exp(-r/a_0) * sin(qx*x[0] + qy*x[1] + qz*x[2]);
}
*/






	
	
