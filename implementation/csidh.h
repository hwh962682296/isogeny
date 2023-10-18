#ifndef CSIDH_H
#define CSIDH_H

#include <stdbool.h>

#include "params.h"

typedef struct private_key {
    int8_t e[NUM_PRIMES]; /* packed int4_t */
} private_key;

typedef struct public_key {
    fp A; /* Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;



extern const public_key base;

void fp_square(fp *x, fp *y,uint64_t* count);
void fp_compute_y(proj *Q,fp *A,fp *y,uint64_t* count);
void montgomery_rhs(fp *rhs, fp const *A, fp const *x,uint64_t* count);/* compute x^3 + Ax^2 + x */
void finding_point(fp *A, proj *P, int i,uint64_t* count);//寻找曲线y^2=x^3+Ax^2+x上阶为l_i的点P，确定型算法
void point_add(proj *P,proj *Q,proj *R,fp *A,uint64_t* count);

void csidh_private(private_key *priv);
bool csidh(public_key *out, public_key const *in, private_key const *priv,uint64_t* count);
void action(public_key *out, public_key const *in, private_key const *priv,uint64_t* count);
void action_one(public_key *out,public_key* in,proj* S,uint64_t* count);

#endif
