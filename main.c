#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

double fun(double x);
double methodDichotomy(double a, double b, double epsilon);

uint16_t main() {
	double x = 1.0;
	double tmp = 0.5 * x*x - 2.0 * x + 2.0;
	printf("%f = %f ?", tmp, fun(x));

	return 0;
}

double fun(const double x) {
	return x < 1.0 ? 1.0 / (sin(x)*sin(x) + 1.0) : 0.5 * x*x - 2.0 * x + 2.0;
}

double methodDichotomy(const double a, const double b, const double epsilon) {
	const double inv2 = 0.5;
	uint16_t k = 0;
	double ak = a;
	double bk = b;
	double ck = (a + b) * inv2;
	const double f_ck = fun(ck);
	
	do {
		double yk = (a + b) * inv2;
		double f_yk = fun(yk);
		
		if(f_yk <= f_ck) {
			bk = ck;
			ck = yk;
		} else {
			double zk = (bk + ck) * inv2;
			double f_zk = fun(zk);
			
			if(f_ck <= f_zk) {
				ak = yk;
				bk = zk;
			} else {
				ak = ck;
				ck = zk;
			}
		}
	} while (fabs(bk - ak) <= epsilon);

	return (ak + bk) * inv2;
}