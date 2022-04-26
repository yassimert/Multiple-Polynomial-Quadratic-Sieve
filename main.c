/*
 * main.c
 *
 *  	Date: June 12, 2017
 *      Author: Mert Yassi
 */
// MULTIPLE POLYNOMIAL QUADRATIC SIEVE
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
#include"mpqs.h"

//mpz_set_str(n, "11145976477", 10);
/*
#define M 1500
#define FBBOUND 600
#define THRESHOLD 5
*/

//mpz_set_str(n, "6241019306901997355512619014111", 10);
#define M 10000
#define FBBOUND 1000
#define THRESHOLD 15

				
//mpz_set_str(n, "26547997890313706964171915281887487041", 10); 
/*#define M 25000
#define FBBOUND 5000
#define THRESHOLD 10
*/

//mpz_set_str(n, "1119348830142138155077694208334174316893", 10);
/*#define M 10000
#define FBBOUND 6000
#define THRESHOLD 10
*/

int main() {
	printf("MULTIPLE POLYNOMIAL QUADRATIC SIEVE\n");
	mpz_t n;
	mpz_init(n);
	mpz_set_str(n, "6241019306901997355512619014111", 10);
	mult_poly_quad_sieve(n, M, FBBOUND, THRESHOLD);
	return 0;
}








