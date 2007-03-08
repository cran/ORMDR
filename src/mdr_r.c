/*
 * mdr_r.c
 *
 * 
 *  by SungGON Yi. <skon@bibs.snu.ac.kr>
 */
#include <R.h>
#include <stdlib.h>
#include <math.h>
/* count the cell number from data */
void make_table(int* c_loc, int* c_num, int lc, 
                int* dataset, int rd, int cd, int* t)
{

	int i, l;
	int loc;
	/* suppose there is no missing data */
	for (i = 0; i < rd; i++) {
		loc = dataset[i];
		for (l = 0; l < lc; l++) 
			loc += c_num[l] * dataset[c_loc[l] + i];
		t[loc]++; 
	}
}
void err_rate(int* comb, int* rcomb, int* ccomb,
              int* train, int* rt, int* ct, int* test, int *rtest,
              double* threshold,
              double* err_train, double* err_test)
{
	int j, k, l;
	/* for calculating error rate */
	int err1, err2;
	int nct = 2 * pow(3, *rcomb);
	double ratio;
	int* c_num = (int*)R_alloc(*rcomb, sizeof(int));
	int* c_loc = (int*) R_alloc(*rcomb, sizeof(int));
	/* counting table */
	int* t_train = (int*) R_alloc(nct, sizeof(int));
	int* t_test = (int*) R_alloc(nct, sizeof(int));
	/* base location of table */
	c_num[0] = 2;
	if (*rcomb > 1)
		for (l = 1; l < *rcomb; l++) 
			c_num[l] = c_num[l - 1] * 3;
	for (j = 0; j < *ccomb; j++) {
		/* initialization */
		for (k = 0; k < nct; k++)
			t_train[k] = t_test[k] = 0;

	/* count the cell from train dataset */ 
		for (l = 0; l < *rcomb; l++)
			c_loc[l] = (comb[*rcomb * j + l] - 1) * *rt;
		make_table(c_loc, c_num, *rcomb, train, *rt, *ct, t_train);
		/* count the cell from test dataset */ 
		for (l = 0; l < *rcomb; l++)
			c_loc[l] = (comb[*rcomb * j + l] - 1) * *rtest;
make_table(c_loc, c_num, *rcomb, test, *rtest, *ct, t_test);
		/* calculate the train and test error rate */
		err1 = err2 = 0;
		for (k = 0; k < nct / 2; k++) {
			/* calculate ratio, declaration is some ambigous */
		if (t_train[k * 2] == 0)
				ratio = *threshold;
			else
				ratio = (double) t_train[k * 2 + 1] / t_train[k * 2]; 
			if (ratio >= *threshold) {
				err1 += t_train[k * 2];
				err2 += t_test[k * 2];
			} else {
				err1 += t_train[k * 2 + 1];
				err2 += t_test[k * 2 + 1];
			}
		}
		err_train[j] = (double) err1 / *rt;
		err_test[j] = (double) err2 / *rtest; 
	}
}

