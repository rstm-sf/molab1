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
double methodDichotomy(double a, double b, double epsilon);
double fun_y(double a, double b);
double methodGoldenSection(double a, double b, double epsilon);
double methodFibonacci(double a, double b, double epsilon);

uint16_t main() {
	double x = 0.5;
	const double epsilon = 1.0e-2;

	//segment_t seg = methodSven(x, 1);
	const double xmin = methodDichotomy(0.5, 3.0, epsilon);

	printf("x = %f\n", xmin);

	system("pause");

	return 0;
}

inline double fun(const double x) {
	return x < 1.0 ? 1.0 / (sin(x)*sin(x) + 1.0) : 0.5 * x*x - 2.0 * x + 2.0;
}

segment_t methodSven(const double x0, const double t) {
	segment_t seg;
	uint16_t k = 0;
	double a = x0 - t;
	double b = x0 + t;
	double x = x0;
	double xk;
	double f_a = fun(a);
	double f_x = fun(x);
	double f_b = fun(b);

	if (f_a <= f_x && f_b <= f_x) {
		printf("Error! No unimodal\n");
		seg.a = 0.0;
		seg.b = 0.0;
		return seg;
	} else if (f_a >= f_x && f_b >= f_x) {
		seg.a = a;
		seg.b = b;
		return seg;
	}

	double delta;

	if (f_a >= f_x && f_x >= f_b ) {
		a = x;
		xk = b;
		delta = t;
	} else if (f_a <= f_x &&  f_x <= f_b  ) {
		b = x;
		xk = a;
		delta = -t;
	}
	k++;

	double f_xk;
	do {
		xk += (2 << k) * delta;
		k++;
		f_xk = fun(xk);

		if (delta == t) { // -->, меняем a; b где=то еще правее
			a = x;
		} else {// <--, меняем b; a где=то еще левее
			b = x;
		}

		x = xk;
	} while (f_xk < f_x);

	if (delta == t) {
		b = xk;
	} else {
		a = xk;
	}

	seg.a = a;
	seg.b = b;

	return seg;
}

double methodDichotomy(const double a, const double b, const double epsilon) {
	const double inv2 = 0.5;
	double ak = a;
	double bk = b;
	double ck = (a + b) * inv2;
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
	return a + (3 - sqrt(5));
}

double methodGoldenSection(const double a, const double b, const double epsilon) {
	double ak = a;
	double bk = b;

	return 0;
}

double methodFibonacci(const double a, const double b, const double epsilon) {
	double ak = a;
	double bk = b;

	return 0;
}