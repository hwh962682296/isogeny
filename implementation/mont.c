
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "uint.h"
#include "fp.h"
#include "mont.h"
#include"gmp.h"

static check(fp* b)
{
    uint c;
    fp_dec(&c,b);
    for(int i=7;i>=0;i--)
    {
        printf("%08X",(uint32_t)(c.c[i]>>32));
        printf("%08X",(uint32_t)(c.c[i]));
    }
    printf("\n\n");
}






static int check_point_eql(proj *P,proj *Q)
{
    fp a,b;
    fp_mul3(&a,&P->x,&Q->z);
    fp_mul3(&b,&P->z,&Q->x);

    return (memcmp(&a,&b,sizeof(fp)));
}







static uint_printf(uint* c)
{
    for(int i=7;i>=0;i--)
    {
        printf("%08X",(uint32_t)(c->c[i]>>32));
        printf("%08X",(uint32_t)(c->c[i]));
    }
    printf("\n\n");
}

//求模逆元 x = a^(-1) mod m
static void extendedGCD(mpz_t a, mpz_t m, mpz_t x, mpz_t y) {
    mpz_t zero, one;
    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(one, 1);

    if (mpz_cmp(m, zero) == 0) {
        mpz_set(x, one);
        mpz_set(y, zero);
        return;
    }

    mpz_t x1, y1;
    mpz_init(x1);
    mpz_init(y1);
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcdext(gcd, x1, y1, m, a);
    
    mpz_set(x, y1);
    mpz_set(y, x1);

    mpz_clear(zero);
    mpz_clear(one);
    mpz_clear(x1);
    mpz_clear(y1);
    mpz_clear(gcd);
}


void CRT(int *crs,unsigned char *message)
{
    mpz_t M, x, Mi, MiInverse, y,result;
    mpz_init(M);
    mpz_init(x);
    mpz_init(Mi);
    mpz_init(MiInverse);
    mpz_init(y);
    mpz_init(result);

    mpz_t a[32],m[32];

    mpz_set_ui(M,1);
    mpz_set_ui(x,0);
    mpz_set_ui(result,0);

    for(int i=0;i<32;i++)
    {
        mpz_init(a[i]);
        mpz_init(m[i]);
        mpz_set_ui(a[i],crs[i]);
        mpz_set_ui(m[i],primes[i+42]);
        mpz_mul(M,M,m[i]);//M=m_1*m_2*...*m_n
    }

    

    for(int i=42;i<NUM_PRIMES;i++)
    {
        mpz_set_ui(Mi,1);
        for(int j=42;j<NUM_PRIMES;j++)
        {
            if(j!=i) mpz_mul(Mi,Mi,m[j-42]);
            else mpz_mul_ui(Mi,Mi,1);
        }//Mi=M/m_i
        extendedGCD(Mi,m[i-42],MiInverse,y);//MiInverse=Mi^(-1) mod m_i
        mpz_mul(x,a[i-42],Mi);
        mpz_mul(x,x,MiInverse);//x=a_i*M_i*t_i
        mpz_add(result,result,x);
        mpz_mod(result,result,M);
    }

    mpz_export(message,NULL,1,1,0,0,result);

    mpz_clear(M);
    mpz_clear(x);
    mpz_clear(Mi);
    mpz_clear(MiInverse);
    mpz_clear(y);
}




//compute P=[m]Q
void Pohlig_Hellman(proj *P, proj *Q, proj *A,uint *m,uint32_t* tag,uint64_t* count)
{
    uint cof;
    uint32_t TAG;
    TAG=tag[0];
    int crs[32]={0};
    for(int i=42;i<NUM_PRIMES;i++)
    {   
        proj P1,Q1;
        uint_set(&cof,1);
        for(int j=42;j<NUM_PRIMES;j++)
        {
            if(j!=i) uint_mul3_64(&cof,&cof,primes[j]);
            else uint_mul3_64(&cof,&cof,1);
        }
        xMUL(&P1,A,P,&cof,count);
        xMUL(&Q1,A,Q,&cof,count);
        if(check_point_eql(&P1,&Q1)==0)
        {
            crs[i-42]=1;
        }
        else
        {
            crs[i-42]=2;
            proj T;
            xDBL(&T,A,&Q1,count);
            proj T1,T2;
            T1.x=Q1.x;
            T1.z=Q1.z; 
            while(check_point_eql(&P1,&T)!=0){
                T2.x=T.x;
                T2.z=T.z;
                xADD(&T,&T,&Q1,&T1,count);
                T1.x=T2.x;
                T1.z=T2.z;
                crs[i-42]++;
            }
        }
        uint32_t mask = 1<<(i-42);
        //printf("%d\n",TAG&mask);
        if((TAG&mask)>0)   crs[i-42]= primes[i]-crs[i-42];
        else crs[i-42]=crs[i-42];
        //printf("%dth crs:%d\n",i-42,crs[i-42]);
    }
    CRT(crs,m);
}



void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A,uint64_t* count)
{
    fp a, b, c, d;

    fp_add3(&a, &Q->x, &Q->z);
    fp_sub3(&b, &Q->x, &Q->z);
    fp_add3(&c, &P->x, &P->z);
    fp_sub3(&d, &P->x, &P->z);
    fp_sq2(&R->x, &c);
    fp_sq2(&S->x, &d);
    fp_mul2(&c, &b);
    fp_mul2(&d, &a);
    fp_sub3(&b, &R->x, &S->x);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_mul3(&R->z, &a, &S->x);
    fp_add3(&S->x, &A->x, &a);
    fp_add2(&R->z, &R->z); /* multiplication by 2 */
    fp_mul2(&R->x, &R->z);
    fp_mul2(&S->x, &b);
    fp_sub3(&S->z, &c, &d);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &c, &d);
    fp_mul2(&R->z, &b);
    fp_sq2(&d, &S->z);
    fp_sq2(&b, &S->x);
    fp_mul3(&S->x, &PQ->z, &b);
    fp_mul3(&S->z, &PQ->x, &d);

    count[0]=count[0]+8;
    count[1]=count[1]+4;
    count[2]=count[2]+11;
}

void xDBL(proj *Q, proj const *A, proj const *P,uint64_t* count)
{
    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);

    count[0]=count[0]+4;
    count[1]=count[1]+2;
    count[2]=count[2]+7;
}

void  xADD(proj *S, proj const *P, proj const *Q, proj const *PQ,uint64_t* count)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);

    count[0]=count[0]+4;
    count[1]=count[1]+2;
    count[2]=count[2]+6;
}

/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k,uint64_t* count)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    unsigned long i = 64 * LIMBS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, A,count);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
int xISOG(proj *A, proj *P, proj const *K, uint64_t k, int check,uint64_t* count)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    fp T[4] = {K->z, K->x, K->x, K->z};
    proj Q;

    fp_mul3(&Q.x,  &P->x, &K->x);
    fp_mul3(&tmp0, &P->z, &K->z);
    fp_sub2(&Q.x,  &tmp0);

    fp_mul3(&Q.z,  &P->x, &K->z);
    fp_mul3(&tmp0, &P->z, &K->x);
    fp_sub2(&Q.z,  &tmp0);
    count[0]=count[0]+4;
    count[2]=count[2]+2;

    proj M[3] = {*K};
    xDBL(&M[1], A, K,count);

    uint64_t i;
    for (i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3],count);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);


        fp_mul3(&tmp0, &P->x, &M[i % 3].x);
        fp_mul3(&tmp1, &P->z, &M[i % 3].z);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.z,  &tmp0);

        count[0]=count[0]+12;
        count[2]=count[2]+4;
    }

    if(check == 1){
        xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3],count);

        proj TestPoint;
        xADD(&TestPoint, &M[i % 3], &M[(i - 1) % 3], K,count);
        
        if (memcmp(&TestPoint.z, &fp_0, sizeof(fp))){
            check = 0;
        } 
    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    count[0] = count[0]+10;
    count[1] = count[1] + 5;
    count[2] = count[2]+5;

    return check;
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
int myxISOG(proj *A, proj *P, int points, proj const *K, uint64_t k, int check,uint64_t* count)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    fp T[4] = {K->z, K->x, K->x, K->z};
    
    proj Q[points];

    for(int p=0; p<points; p++){
        fp_mul3(&Q[p].x,  &P[p].x, &K->x);
        fp_mul3(&tmp0, &P[p].z, &K->z);
        fp_sub2(&Q[p].x,  &tmp0);

        fp_mul3(&Q[p].z,  &P[p].x, &K->z);
        fp_mul3(&tmp0, &P[p].z, &K->x);
        fp_sub2(&Q[p].z,  &tmp0);
        count[0]=count[0]+4;
        count[2]=count[2]+2;
    }

    proj M[3] = {*K};
    xDBL(&M[1], A, K,count);

    uint64_t i;
    for (i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3],count);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);
        count[0]=count[0]+6;
        count[2]=count[2]+2;

        for(int p=0; p<points ; p++){
            fp_mul3(&tmp0, &P[p].x, &M[i % 3].x);
            fp_mul3(&tmp1, &P[p].z, &M[i % 3].z);
            fp_sub2(&tmp0, &tmp1);
            fp_mul2(&Q[p].x,  &tmp0);

            fp_mul3(&tmp0, &P[p].x, &M[i % 3].z);
            fp_mul3(&tmp1, &P[p].z, &M[i % 3].x);
            fp_sub2(&tmp0, &tmp1);
            fp_mul2(&Q[p].z,  &tmp0);
            count[0]=count[0]+6;
            count[2]=count[2]+2;
        }
    }

    if(check == 1){
        xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3],count);

        proj TestPoint;
        xADD(&TestPoint, &M[i % 3], &M[(i - 1) % 3], K,count);
        
        if (memcmp(&TestPoint.z, &fp_0, sizeof(fp))){
            check = 0;
        } 
    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    count[0]=count[0]+8;
    count[1]=count[1]+3;
    count[2]=count[2]+5;

    /* X := X * Xim^2, Z := Z * Zim^2 */
    for (int i = 0; i < points; ++i)
    {
        fp_sq1(&Q[i].x);
        fp_sq1(&Q[i].z);
        fp_mul2(&P[i].x, &Q[i].x);
        fp_mul2(&P[i].z, &Q[i].z);

        count[0]=count[0]+2;
        count[1]=count[1]+2;
    }

    return check;
}

