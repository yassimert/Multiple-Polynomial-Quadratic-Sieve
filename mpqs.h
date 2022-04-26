/*
 * mpqs.h
 *
 *  	Date: June 12, 2017
 *      Author: Mert Yassi
 */
// MULTIPLE POLYNOMIAL QUADRATIC SIEVE

#ifndef MPQS_H_
#define MPQS_H_

typedef signed long long int ull_int;

typedef struct ARRAYS *ARRAYS_PT;
typedef struct ARRAYS{
	ull_int numOfPrimes;
	ull_int *factors, *exponents;
}ARRAYS_T[1];

static ull_int ctrSize = 0;

void sieve_of_eratosthenes(ARRAYS_PT ARRS, int n);
void mult_poly_quad_sieve(mpz_t n, ull_int M, ull_int FBBOUND, ull_int THRESHOLD);
void calc_nullspace(ARRAYS_PT ARRS, ull_int **arr, ull_int fbsize);
void evaluate_nullspace(mpz_t n, mpz_t *kArr, mpz_t *lArr);
int mpz_sqrt_mod_p(mpz_t q, const mpz_t n, const mpz_t p);
int mpz_sqrt_mod_pe(mpz_t q, const mpz_t a, const mpz_t p, const int e);

#endif /* MPQS_H_ */

