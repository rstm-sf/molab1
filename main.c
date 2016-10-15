#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

struct segment {
	double a;
	double b;
};

double fun(double x);
struct segment methodSven(double x0, double t);
double methodDichotomy(double a, double b, double epsilon);

uint16_t main() {
	double x = 1.0;

	struct segment seg = methodSven(x, 0.1);
	x = methodDichotomy(seg.a, seg.b, 1.0e-2);

	printf("x = %f\n", x);

	return 0;
}

double fun(const double x) {
	return x < 1.0 ? 1.0 / (sin(x)*sin(x) + 1.0) : 0.5 * x*x - 2.0 * x + 2.0;
}

struct segment methodSven(const double x0, const double t) {
	struct segment seg;
	uint16_t k = 0;
	double a = x0 - t;
	double b = x0 + t;
	double x = x0;
	double xk;
	double f_xmt = fun(a);
	double f_x = fun(x);
	double f_xpt = fun(b);

	if (f_xmt <= f_x && f_xpt <= f_x) {
		printf("Error! No unimodal\n");
		seg.a = 0.0;
		seg.b = 0.0;
		return seg;
	} else if (f_xmt >= f_x && f_xpt >= f_x) {
		seg.a = a;
		seg.b = b;
		return seg;
	}

	double delta;

	if (f_xmt >= f_x && f_xpt <= f_x) {
		a = x;
		xk = x + t;
		delta = t;
	} else if (f_xmt <= f_x && f_xpt >= f_x) {
		b = x;
		xk = x - t;
		delta = -t;
	}
	k++;

	double f_xk;
	do {
		xk += (2 << k) * delta;
		k++;
		f_xk = fun(xk);

		if (delta == t) {
			a = x;
		} else if (delta == -t) {
			b = x;
		}
	} while (f_xk < f_x);


	if (delta == t) {
		b = xk;
	} else if (delta == -t) {
		a = xk;
	}

	seg.a = a;
	seg.b = b;

	return seg;
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

		if (f_yk <= f_ck) {
			bk = ck;
			ck = yk;
		}
		else {
			double zk = (bk + ck) * inv2;
			double f_zk = fun(zk);

			if (f_ck <= f_zk) {
				ak = yk;
				bk = zk;
			}
			else {
				ak = ck;
				ck = zk;
			}
		}
	} while (fabs(bk - ak) <= epsilon);

	return (ak + bk) * inv2;
}