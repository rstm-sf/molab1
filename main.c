#include <stdio.h>
#include <stdlib.h>
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
double methodGoldenSection(segment_t seg, double epsilon);
uint32_t findNforFibonacci(segment_t seg, double l);
double methodFibonacci(segment_t seg, double l, double epsilon);

uint32_t main() {
	const double x0 = 1.5;
	const double epsilon = 1.0e-2;
	const double t = 0.1;
	const double l = epsilon;

	const segment_t seg = methodSven(x0, t);
	printf("[%.2f, %.2f]\n", seg.a, seg.b);
	const double xmin = methodFibonacci(seg, l, epsilon);

	printf("\nxmin = %f\n", xmin);

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
		k++;
		xk = xk1;
		f_xk = fun(xk);
		xk1 = xk + (1 << k) * delta;
		f_xk1 = fun(xk1);
	} while (f_xk1 < f_xk);

	*bound_xk = xk;
	*bound_xk1 = xk1;

	seg.a = a;
	seg.b = b;

	printf("The method Sven\nIters = %" PRIu16 "\n", k);

	return seg;
}

double methodDichotomy(const segment_t seg, const double epsilon) {
	const double inv2 = 0.5;
	double ak = seg.a;
	double bk = seg.b;
	double ck = (ak + bk) * inv2;
	double f_ck = fun(ck);
	double convergence = fabs(bk - ak);

	uint16_t k = 0;
	while (convergence > epsilon) {
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
		convergence = fabs(bk - ak);
		k++;
	}

	printf("The method of Dichotomy\nIters = %" PRIu16 "\nConvergence = %e\n", k, convergence);

	return (ak + bk) * inv2;
}

double methodGoldenSection(const segment_t seg, const double epsilon) {
	double ak = seg.a;
	double bk = seg.b;
	// (3 - sqrt(5))/2 ~ 0.3819660112501051
	double yk = ak + 0.3819660112501051 * (bk - ak);
	double zk = ak + bk - yk;
	double f_yk = fun(yk);
	double f_zk = fun(zk);
	double convergence;

	uint16_t k = 0;
	do {
		if (f_yk <= f_zk) {
			bk = zk;
			zk = yk;
			f_zk = f_yk;
			yk = ak + bk - yk;
			f_yk = fun(yk);
		} else {
			ak = yk;
			yk = zk;
			f_yk = f_zk;
			zk = ak + bk - zk;
			f_zk = fun(zk);
		}

		convergence = fabs(ak - bk);
		k++;
	} while (convergence > epsilon);

	printf("The method of the Golden section\nIters = %" PRIu16 "\nConvergence = %e\n", k, convergence);

	return (ak + bk) * 0.5;
}

uint32_t findNforFibonacci(const segment_t seg, const double l) {
	uint32_t N = 1;
	uint64_t F0 = 1, F1 = 2, F2 = 3;
	const uint64_t bound = (uint64_t)(fabs(seg.b - seg.a) / l);

	while (F2 < bound) {
		F0 = F1;
		F1 = F2;
		F2 = F1 + F0;
		N++;
	}

	return N;
}

double methodFibonacci(const segment_t seg, const double l, const double epsilon) {
	double ak = seg.a;
	double bk = seg.b;

	const uint32_t N = findNforFibonacci(seg, l); // N >= 1
	double *F = (double *)malloc((N + 3) * sizeof(double));
	F[0] = F[1] = 1.0; F[2] = 2.0; F[3] = 3.0;
	double *currentF = F + 1;
	for (uint32_t i = 4; i < N + 3; ++i) {
		currentF++;
		currentF[2] = currentF[1] + currentF[0];
	}

	double yk = ak + currentF[0] / currentF[2] * (bk - ak);
	double zk = ak + bk - yk;
	double f_yk = fun(yk);
	double f_zk = fun(zk);

	uint32_t k = 0;
	while (k < N) {
		k++;
		currentF--;
		if (f_yk <= f_zk) {
			bk = zk;
			zk = yk;
			f_zk = f_yk;
			yk = ak + currentF[0] / currentF[2] * (bk - ak);
			f_yk = fun(yk);
		} else {
			ak = yk;
			yk = zk;
			f_yk = f_zk;
			zk = ak + currentF[1] / currentF[2] * (bk - ak);
			f_zk = fun(zk);
		}
	}

	zk += epsilon;
	f_zk = fun(zk);
	if (f_yk <= f_zk) {
		bk = zk;
	}

	printf("The method Fibonacci\nN = %" PRIu32 "\n", N);
	free(F);

	return (ak + bk) * 0.5;
}