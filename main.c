#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

double fun(double x);

uint16_t main() {
	double x = 1.0;
	double tmp = 0.5 * x*x - 2.0 * x + 2.0;
	printf("%f = %f ?", tmp, fun(x));

	return 0;
}

double fun(const double x) {
	return x < 1.0 ? 1.0 / (sin(x)*sin(x) + 1.0) : 0.5 * x*x - 2.0 * x + 2.0;
}