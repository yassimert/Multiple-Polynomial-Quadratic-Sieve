/*
 * mpqs.c
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

#define SIZE 1000000

void sieve_of_eratosthenes(ARRAYS_PT ARRS, int n) {
	ull_int *primes, *factorArr, i, j, count = 0, k = 0;
	primes = (ull_int*)malloc(sizeof(ull_int) * n);

	for(i = 2; i < n; i++){
		primes[i] = 1; // make all elements true
	}

	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			for(j = i; (i * j) < n; j++){
				primes[i * j] = 0; // mark composites
			}
		}
	}

	for(i = 2; i < n; i++){ // loop for determining numOfPrimes
		if(primes[i] == 1){
			count++;
		}
	}

	factorArr = (ull_int*)malloc(sizeof(ull_int) * count);
	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			factorArr[k] = i;
			k++;
		}
	}
	ARRS->factors = factorArr;
	ARRS->numOfPrimes = count;
	free(primes);
}

void mult_poly_quad_sieve(mpz_t n, ull_int M, ull_int FBBOUND, ull_int THRESHOLD) {
	ARRAYS_PT ARRS = (ARRAYS_PT)malloc(sizeof(ARRAYS_T));
	sieve_of_eratosthenes(ARRS, FBBOUND);
	//******** INITIALIZE START ********//
	ull_int i, j, r, k, aj = 0, fbsize = 0, control2 = 0, SIb, SIfbase, lSIfbase, ctr = 0, tsb, ttsb, SIainv, mainv, mmainv, soln1, soln2;
	ull_int *indicesArr, **arr, *tLog, *tSqrt, *sieveArr;
	mpz_t mSqrt, q, mn, rn, prime, x1, a, b, ainv, ax, axpb, axpbs, pol, kTemp, modTemp;
	mpz_t *mpzfactors, *mpzfbase, *axbArr, *axbsArr, *pols, *uArr, *vArr;
	ctrSize = 0;

	mpz_init(mSqrt);
	mpz_init(q);
	mpz_init(mn);
	mpz_init(rn);
	mpz_init(prime);
	mpz_init(x1);
	mpz_init(a);
	mpz_init(b);
	mpz_init(ainv);
	mpz_init(ax);
	mpz_init(axpb);
	mpz_init(axpbs);
	mpz_init(pol);
	mpz_init(kTemp);
	mpz_init(modTemp);

	mpzfactors = (mpz_t*)malloc(sizeof(mpz_t) * ARRS->numOfPrimes);
	for(i = 0; i < ARRS->numOfPrimes; i++){
		mpz_init_set_si(mpzfactors[i], ARRS->factors[i]);
	}

	sieveArr = (ull_int*)malloc(sizeof(ull_int) * M);
	ull_int *pt = &sieveArr[M/2];

	for(i = -M/2; i < M/2; i++){
		pt[i] = 0;
	}

	for(i = 0; i < ARRS->numOfPrimes; i++){
		if(mpz_jacobi(n, mpzfactors[i]) == 1){ // create bound for factor base initialization
			fbsize++;
		}
	}

	vArr = (mpz_t*)malloc(sizeof(mpz_t) * SIZE);
	uArr = (mpz_t*)malloc(sizeof(mpz_t) * SIZE);
	for(i = 0; i < SIZE; i++){
		mpz_init(vArr[i]);
		mpz_init(uArr[i]);
	}

	//******** INITIALIZE END ********//

	j = 0;
	if(mpz_jacobi(n, mpzfactors[0]) == 1){ // determine the factor base primes p < f such that the jacobi symbol (n/p) is 1
		mpzfbase = (mpz_t*)malloc(sizeof(mpz_t) * fbsize);
		for(i = 0; i < fbsize; i++){
			mpz_init(mpzfbase[i]);
		}
		for(i = 0; i < ARRS->numOfPrimes; i++){
			if(mpz_jacobi(n, mpzfactors[i]) == 1){
				mpz_set(mpzfbase[j], mpzfactors[i]); // choose quadratic residues
				j++;
			}
		}
	}else{
		j++;
		fbsize = fbsize + 1;
		mpzfbase = (mpz_t*)malloc(sizeof(mpz_t) * fbsize);
		for(i = 0; i < fbsize; i++){
			mpz_init(mpzfbase[i]);
		}
		mpz_set(mpzfbase[0], mpzfactors[0]); // add 2 to factor base
		for(i = 1; i < ARRS->numOfPrimes; i++){
			if(mpz_jacobi(n, mpzfactors[i]) == 1){
				mpz_set(mpzfbase[j], mpzfactors[i]); // choose quadratic residues
				j++;
			}
		}

	}
	printf("Factorbase: ");
	for(i = 0; i < fbsize; i++){
		mpz_out_str(stdout, 10, mpzfbase[i]);
		printf(" ");
	}
	printf("\n");
	arr = (ull_int**)malloc(sizeof(ull_int*) * fbsize);
	for(i = 0; i < fbsize; i++){
		arr[i] = (ull_int*)malloc(sizeof(ull_int) * SIZE);
	}
	tLog = (ull_int*)malloc(sizeof(ull_int) * fbsize);
	tSqrt = (ull_int*)malloc(sizeof(ull_int) * fbsize);
	ARRS->exponents = (ull_int*)malloc(sizeof(ull_int) * fbsize);

	for(i = 0; i < fbsize; i++){ // for each factor base prime p,
		mpz_sqrt_mod_p(mSqrt, n, mpzfbase[i]); // compute mSqrt, a modular square root of n (mod p)
		tSqrt[i] = mpz_get_si(mSqrt); // store mSqrt
	}

	for(i = 0; i < fbsize; i++){
		SIfbase = mpz_get_si(mpzfbase[i]);
		lSIfbase = round(log2(SIfbase)); //compute the base-2 logarithm of each prime lSIfbase = ⌊log2 p⌉(rounding)
		tLog[i] = lSIfbase;	// store ⌊log2 p⌉
	}

	// find a prime q ≈ √( √(2n) / m) such that the jacobi symbol (n/q) is 1
	mpz_mul_si(mn, n, 2); // 2*n
	mpz_sqrt(rn, mn); // √(2*n)
	mpz_cdiv_q_ui(rn, rn, M/2); // √(2*n)/(M/2)
	mpz_sqrt(q, rn); // √(√(2*n)/(M/2))

	mpz_nextprime(prime, q); // check if n is quadratic residue mod prime otherwise choose another prime
	while(mpz_jacobi(n, prime) != 1){
		mpz_nextprime(prime, prime);
	}
	while(control2 != 1){
		for(i = -M/2; i < M/2; i++){
			pt[i] = 0;
		}
		mpz_sqrt_mod_p(x1, n, prime); // x1^2 = n (mod prime) find modular square root
		mpz_mul(a, prime, prime); // a = prime^2
		mpz_sqrt_mod_pe(b, n, prime, 2); // b^2 = n (mod prime^2) find modular square root via hensel lifting
		SIb = mpz_get_si(b);

		//******** compute solutions START ********//
		for(i = 1; i < fbsize; i++){ // for each odd prime p in the factor base, compute soln1,soln2 = a^-1((-)tSqrt[p] − b) % mpzfbase[p]
			mpz_invert(ainv, a, mpzfbase[i]); // a^-1
			SIainv = mpz_get_si(ainv);
			tsb = tSqrt[i] - SIb; // (tSqrt[i] - b)
			mainv = SIainv * tsb; // (a^-1 * (tSqrt[i] - b))
			SIfbase = mpz_get_si(mpzfbase[i]);
			soln1 = mainv % SIfbase; // soln1 =  (a^-1 * (tSqrt[i] - b)) % mpzfbase[i]
			j = 0;
			while((soln1 + (j * SIfbase) > -M/2) && (soln1 + (j * SIfbase) < M/2)){
				pt[soln1 + (j * SIfbase)] = pt[soln1 + (j * SIfbase)] + tLog[i];
				j++;
			}
			ttsb = -tSqrt[i] - SIb; // (-tSqrt[i] - b)
			mmainv = SIainv * ttsb; // (a^-1 * (-tSqrt[i] - b))
			soln2 =  mmainv % SIfbase; // soln2 =  (a^-1 * (-tSqrt[i] - b)) % mpzfbase[i]
			j = 0;
			while((soln2 + (j * SIfbase) > -M/2) && (soln2 + (j * SIfbase) < M/2)){
				pt[soln2 + (j * SIfbase)] = pt[soln2 + (j * SIfbase)] + tLog[i];
				j++;
			}
			//******** compute solutions END ********//

			ctr = 0;
			for(i = -M/2; i < M/2; i++){
				if(pt[i] > THRESHOLD){ // > log(m√n) - small error
					ctr++;
				}
			}
			indicesArr = (ull_int*)malloc(sizeof(ull_int) * ctr);
			j = 0;
			for(i = -M/2; i < M/2; i++){
				if(pt[i] > THRESHOLD){ // > log(m√n) - small error
					indicesArr[j] = i;
					j++;
				}
			}

			//******** create polynomials START ********//
			axbArr = (mpz_t*)malloc(sizeof(mpz_t) * ctr);
			axbsArr = (mpz_t*)malloc(sizeof(mpz_t) * ctr);
			pols = (mpz_t*)malloc(sizeof(mpz_t) * ctr);
			for(i = 0; i < ctr; i++){
				mpz_init(axbArr[i]);
				mpz_init(axbsArr[i]);
				mpz_init(pols[i]);
			}
			for(i = 0; i < ctr; i++){
				mpz_mul_si(ax, a, indicesArr[i]); // a*x
				mpz_add(axpb, ax, b);  // a*x + b
				mpz_set(axbArr[i], axpb);
				mpz_mul(axpbs, axpb, axpb);  // (a*x + b)^2
				mpz_sub(pol, axpbs, n); //  (a*x + b)^2 - n
				mpz_set(axbsArr[i], pol);
				mpz_divexact(pols[i], pol, a); //  ((a*x + b)^2 - n) / a
			}
			//******** create polynomials END ********//
			k = 0;
			for(r = 0; r < ctr; r++){
				for(i = 0; i < fbsize; i++){
					ARRS->exponents[i] = 0;
				}
				int control = 0; i = 0;
				while(control == 0){
					mpz_mod(modTemp, pols[k], mpzfbase[i]);
					if((mpz_sgn(modTemp) == 0) && (mpz_sgn(pols[k]) != 0)){
						while((mpz_sgn(modTemp) == 0) && (mpz_sgn(pols[k]) != 0)){
							mpz_divexact(pols[k], pols[k], mpzfbase[i]);
							ARRS->exponents[i]++;
							mpz_mod(modTemp, pols[k], mpzfbase[i]);
						}
					}else{
						if(i < fbsize - 1){
							i++;
						}else{
							control = 1;
						}
					}
				}
				if(mpz_cmp_ui(pols[k], 1) == 0){
					mpz_set(vArr[aj], axbsArr[k]);
					mpz_set(uArr[aj], axbArr[k]);
					printf(" [ ");
					for(i = 0; i < fbsize; i++){
						arr[i][aj] = ARRS->exponents[i] % 2;
						printf("%lld ",arr[i][aj]);
					}
					printf(" ]\n");
					ctrSize++;
					aj++;
				}
				k++;
			}

			if(ctrSize > ((fbsize*105)/100)){
				control2 = 1;
			}else{
				mpz_nextprime(prime, prime);
				while(mpz_jacobi(n, prime) != 1){
					mpz_nextprime(prime, prime);
				}
			}
			free(indicesArr);
			for(i = 0; i < ctr; i++){
				mpz_clear(pols[i]);
				mpz_clear(axbArr[i]);
				mpz_clear(axbsArr[i]);
			}
		}
	}

	calc_nullspace(ARRS, arr, fbsize);
	evaluate_nullspace(n, uArr, vArr);

	//******** FREE START ********//

	for(i = 0; i < ARRS->numOfPrimes; i++){
		mpz_clear(mpzfactors[i]);
	}

	for(i = 0; i < fbsize; i++){
		mpz_clear(mpzfbase[i]);
	}

	for(i = 0; i < SIZE; i++){
		mpz_clear(vArr[i]);
		mpz_clear(uArr[i]);
	}

	mpz_clear(mSqrt);
	mpz_clear(q);
	mpz_clear(mn);
	mpz_clear(rn);
	mpz_clear(prime);
	mpz_clear(x1);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(ainv);
	mpz_clear(ax);
	mpz_clear(axpb);
	mpz_clear(axpbs);
	mpz_clear(pol);
	mpz_clear(kTemp);
	mpz_clear(modTemp);

	free(tLog);
	free(tSqrt);
	free(sieveArr);
	free(ARRS->factors);
	free(ARRS->exponents);
	for(i = 0; i < fbsize; i++){
		free(arr[i]);
	}
	free(arr);
	free(ARRS);
	//******** FREE END ********//
}

void calc_nullspace(ARRAYS_PT ARRS, ull_int **arr, ull_int fbsize) {
	ull_int i, j, count = 0, count2 = 0;
	FILE *f = fopen("matrix.txt", "w");
	if (f == NULL){
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f,"M:=Matrix(GF(2),%lld,%lld,[", ctrSize, fbsize);
	for(j = 0; j < ctrSize; j++){
		count++;
		for(i = 0; i < fbsize; i++){
			count2++;
			fprintf(f,"%lld",arr[i][j]);
			if(count2 < fbsize){
				fprintf(f,",");
			}
		}
		if(count < ctrSize){
			fprintf(f,",");
		}
		count2 = 0;
	}
	fprintf(f,"]);");
	fprintf(f,"\nN:=NullSpace(M);");
	fprintf(f,"\nB:=Basis(N);");
	fprintf(f,"\nB;");
	fclose(f);
	//int ret = system("~/magma < matrix.txt > result.txt"); // you can factor larger numbers using local Magma
	int ret = system("python3 encoder.py matrix.txt > result.txt"); // will not work if the size of matrix.txt > 50 kB (Online Magma input limit)
}

void evaluate_nullspace(mpz_t n, mpz_t *kArr, mpz_t *lArr) {
	mpz_t old_x, old_y, x, y2, y, xpy, xmy, p, q;
	ull_int i, j, k;
	mpz_init(old_x);
	mpz_init(old_y);
	mpz_init(y);
	mpz_init(xpy);
	mpz_init(xmy);
	mpz_init(p);
	mpz_init(q);
	mpz_init_set_ui(x, 1);
	mpz_init_set_ui(y2, 1);

	FILE *f = fopen("result.txt", "rb");
	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);
	char *string = malloc(fsize + 1);
	int ret = fread(string, fsize, 1, f);
	fclose(f);
	string[fsize] = 0;

	int sct = 1; // sct = 0 when using local Magma
	i = 0;
	while(string[i] != '\0'){
		if(string[i] == 44){
			sct++;
		}
		i++;
	}

	char **arr2;
	arr2 = (char**)malloc(sizeof(char*) * sct);
	for(i = 0; i < sct; i++){
		arr2[i] = (char*)malloc(sizeof(char) * ctrSize);
	}

	i = 0, j = 0, k = 0;
	while(string[i] != '\0'){
		if(string[i] == '('){
			while(string[i] != ')'){
				if(((string[i] == 48) || (string[i] == 49))){
					arr2[k][j] = string[i];
					j++;
				}
				i++;
			}
			j = 0;
			k++;
		}else{
			i++;
		}
	}

	if(sct == 1){
		printf("Magma could not find NullSpace of M!\n");
	}else{
		printf("NullSpace:\n");
		for(i = 0; i < sct; i++){
			printf(" [ ");
			for(j = 0; j < ctrSize; j++){
				printf("%c ",arr2[i][j]);
			}
			printf("]\n");
		}
	}
	int ctrl = 0;
	for(i = 0; i < sct; i++){
		for(j = 0; j < ctrSize; j++){
			if(arr2[i][j] == 49){
				mpz_set(old_x, kArr[j]);
				mpz_set(old_y, lArr[j]);
				mpz_mul(x, old_x, x);
				mpz_mul(y2, old_y, y2);
			}
		}
		mpz_sqrt(y, y2);
		mpz_add(xpy, x, y);
		mpz_sub(xmy, x, y);
		mpz_gcd(p, xpy, n);
		mpz_gcd(q, xmy, n);

		if(((mpz_cmp_ui(p, 1) != 0) && (mpz_cmp(p, n) != 0)) && ((mpz_cmp_ui(q, 1) != 0) && (mpz_cmp(q, n) != 0))){
			printf("1ST FACTOR P: ");
			mpz_out_str(stdout, 10, p);
			printf("\n");
			printf("2ND FACTOR Q: ");
			mpz_out_str(stdout, 10, q);
			printf("\n");
			ctrl = 1;
			i = sct;
		}else{
			mpz_set_ui(x, 1);
			mpz_set_ui(y2, 1);
		}
	}
	if(ctrl != 1){
		printf("The algorithm could not find proper factors!\n");
	}

}

int mpz_sqrt_mod_p(mpz_t q, const mpz_t n, const mpz_t p) {
	mpz_t w, n_inv, y;
	unsigned int i, s;

	if(mpz_divisible_p(n, p)){ // Is n a multiple of p?            
		mpz_set_ui(q, 0); // Yes, then the square root is 0.  
		return 1; // (special case, not caught        
	} // otherwise)                       
	if(mpz_legendre(n, p) != 1) // Not a quadratic residue?         
		return 0; // No, so return error              
	if(mpz_tstbit(p, 1) == 1){ // p = 3 (mod 4) ?                  
		mpz_set(q, p);
		mpz_add_ui(q, q, 1);
		mpz_fdiv_q_2exp(q, q, 2);
		mpz_powm(q, n, q, p); // q = n ^ ((p+1) / 4) (mod p)      
		return 1;
	}
	mpz_init(y);
	mpz_init(w);
	mpz_init(n_inv);
	mpz_set(q, p);
	mpz_sub_ui(q, q, 1); // q = p-1                          
	s = 0; // Factor out 2^s from q            
	while(mpz_tstbit(q, s) == 0)
		s++;
	mpz_fdiv_q_2exp(q, q, s); // q = q / 2^s                      
	mpz_set_ui(w, 2); // Search for a non-residue mod p   
	while(mpz_legendre(w, p) != -1) // by picking the first w such that 
		mpz_add_ui(w, w, 1); // (w/p) is -1                      
	mpz_powm(w, w, q, p); // w = w^q (mod p)                  
	mpz_add_ui(q, q, 1);
	mpz_fdiv_q_2exp(q, q, 1); // q = (q+1) / 2                    
	mpz_powm(q, n, q, p); // q = n^q (mod p)                  
	mpz_invert(n_inv, n, p);
	for(;;){
		mpz_powm_ui(y, q, 2, p); // y = q^2 (mod p)                  
		mpz_mul(y, y, n_inv);
		mpz_mod(y, y, p); // y = y * n^-1 (mod p)             
		i = 0;
		while(mpz_cmp_ui(y, 1) != 0){
			i++;
			mpz_powm_ui(y, y, 2, p); //  y = y ^ 2 (mod p)
			 
		}
		if(i == 0){ // q^2 * n^-1 = 1 (mod p), return   
			mpz_clear(n_inv);
			mpz_clear(w);
			mpz_clear(y);
			return 1;
		}
		if(s - i == 1){ // In case the exponent to w is 1, 
			mpz_mul(q, q, w); // Don't bother exponentiating     
		}else{
			mpz_powm_ui(y, w, 1 << (s - i - 1), p);
			mpz_mul(q, q, y);
		}
		mpz_mod(q, q, p); // r = r * w^(2^(s-i-1)) (mod p)  
	}
	return 0;
}

int mpz_sqrt_mod_pe(mpz_t q, const mpz_t a, const mpz_t p, const int e) {
	int r, i;
	mpz_t abb, b, t, pf, pf1, b2, h;

	mpz_init(h);
	mpz_init(b2);
	mpz_init(pf);
	mpz_init(pf1);
	mpz_init(b);
	mpz_init(t);
	mpz_init(abb);

	r = mpz_sqrt_mod_p(b, a, p);

	// Hensel lifting
	mpz_set(pf, p);
	for(i = 1; i < e; i++){
		mpz_mul(pf1, pf, p);
		mpz_mul(abb, b, b); // b^2
		mpz_sub(abb, a, abb); // a-b^2
		mpz_divexact(abb, abb, pf); // (a-b^2) div p^f 
		mpz_mul_ui(t, b, 2); // 2*b
		mpz_invert(t, t, pf1); // (2*b)^(-1) mod p^(f+1) 
		mpz_mul(h, abb, t); // (2*b)^(-1) * [ (a-b^2) div p^f ] 
		mpz_mod(h, h, pf1); // [ (2*b)^(-1) * [ (a-b^2) div p^f ] ] mod p^(f+1) 
		mpz_mul(h, h, pf);
		mpz_add(b, b, h);
		mpz_mod(b, b, pf1);
		mpz_set(pf, pf1);
	}
	mpz_set(q, b);

	mpz_clear(h);
	mpz_clear(b2);
	mpz_clear(pf);
	mpz_clear(pf1);
	mpz_clear(b);
	mpz_clear(t);
	mpz_clear(abb);
	return r;
}


