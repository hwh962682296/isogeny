#ifndef MONT_H
#define MONT_H

#include "params.h"

void CRT(int* crs,unsigned char* message);
void Pohlig_Hellman(proj *P, proj *Q, proj *A,uint* m,uint32_t* tag,uint64_t* count);
void xDBL(proj *Q, proj const *A, proj const *P, uint64_t* count);
void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ,uint64_t* count);
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A,uint64_t* count);
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k,uint64_t* count);
int xISOG(proj *A, proj *P, proj const *K, uint64_t k, int check,uint64_t* count);
int myxISOG(proj *A, proj *P, int points, proj const *K, uint64_t k, int check,uint64_t* count);

#endif
