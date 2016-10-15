#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

struct segment {
	double a;
	double b;
};

typedef struct segment segment_t;

double fun(double x);
segment_t methodSven(double x0, double t);
double methodDichotomy(segment_t seg, double epsilon);
double fun_y(double a, double b);
double methodGoldenSection(segment_t seg, double epsilon);
double methodFibonacci(segment_t seg, double epsilon);

uint16_t main() {
	double x0 = 3.0;
	const double epsilon = 1.0e-2;
	const double t = 0.1;

	segment_t seg = methodSven(x0, t);
	printf("[%.2f, %.2f]\n", seg.a, seg.b);
	const double xmin = methodDichotomy(seg, epsilon);

	printf("x = %f\n", xmin);

	system("pause");

	return 0;
}

inline double fun(const double x) {
	return x < 1.0 ? 1.0 / (sin(x*x) + 1.0) : 0.5 * x*x - 2.0 * x + 2.0;
}

segment_t methodSven(const double x0, const double t) {
	segment_t seg = { 0 };
	uint16_t k = 0;
	double a = x0 - t;
	double b = x0 + t;
	const double f_a = fun(a);
	double f_xk = fun(x0);
	const double f_b = fun(b);

	if (f_a <= f_xk && f_b <= f_xk) {
		printf("Error! No unimodal\n");
		return seg;
	} else if (f_a >= f_xk && f_b >= f_xk) {
		seg.a = a;
		seg.b = b;
		return seg;
	}

	double delta, xk1;
	double *bound_xk, *bound_xk1;
	if (f_a >= f_xk && f_xk >= f_b) {
		a = x0;
		xk1 = b;
		delta = t;
		bound_xk = &a;
		bound_xk1 = &b;
	} else {
		b = x0;
		xk1 = a;
		delta = -t;
		bound_xk1 = &a;
		bound_xk = &b;
	}

	double f_xk1, xk;
	do {
		xk = xk1;
		f_xk = fun(xk);
		xk1 = xk + (2 << k) * delta;
		f_xk1 = fun(xk1);
		k++;
	} while (f_xk1 < f_xk);

	*bound_xk = xk;
	*bound_xk1 = xk1;

	seg.a = a;
	seg.b = b;

	return seg;
}

double methodDichotomy(const segment_t seg, const double epsilon) {
	const double inv2 = 0.5;
	double ak = seg.a;
	double bk = seg.b;
	double ck = (ak + bk) * inv2;
	double f_ck = fun(ck);
	double convergense = fabs(bk - ak);

	while (convergense > epsilon) {
		const double yk = (ak + bk) * inv2;
		const double f_yk = fun(yk);

		if (f_yk <= f_ck) {
			bk = ck;
			ck = yk;
			f_ck = fun(ck);
		} else {
			const double zk = (bk + ck) * inv2;
			const double f_zk = fun(zk);

			if (f_ck <= f_zk) {
				ak = yk;
				bk = zk;
			} else {
				ak = ck;
				ck = zk;
				f_ck = fun(ck);
			}
		}
		convergense = fabs(bk - ak);
	}

	return (ak + bk) * inv2;
}

inline double fun_y(const double a, const double b) {
	// sqrt(5) ~ 2.2360679774997898
	return a + (3 - 2.2360679774997898);
}

double methodGoldenSection(const segment_t seg, const double epsilon) {
	double ak = seg.a;
	double bk = seg.b;

	return 0;
}

double methodFibonacci(const segment_t seg, const double epsilon) {
	double ak = seg.a;
	double bk = seg.b;

	return 0;
}