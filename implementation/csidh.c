
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "uint.h"
#include "fp.h"
#include "mont.h"
#include "csidh.h"
#include "rng.h"
#include "classgroup.h"

//#define ORIGINAL
#define UNIFORM

const public_key base = {{{0}}}; /* A = 0 */

#ifdef UNIFORM
void csidh_private(private_key *priv)
{
    sample_from_classgroup(priv->e);
}
#endif

#ifdef ORIGINAL
#define MAX_EXPONENT 5
/* TODO allow different encodings depending on parameters */
/* TODO waste less randomness */
void csidh_private(private_key *priv)
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < NUM_PRIMES; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->e[i] = buf[j];
                if (++i >= NUM_PRIMES)
                    break;
            }
        }
    }
}
#endif

static bool validate_rec(proj *P, proj const *A, size_t lower, size_t upper, uint *order, bool *is_supersingular,uint64_t* count)
{
    assert(lower < upper);

    if (upper - lower == 1) {

        /* now P is [(p+1) / l_lower] times the original random point */
        /* we only gain information if this multiple is non-zero */

        if (memcmp(&P->z, &fp_0, sizeof(fp))) {

            uint tmp;
            uint_set(&tmp, primes[lower]);
            xMUL(P, A, P, &tmp,count);

            if (memcmp(&P->z, &fp_0, sizeof(fp))) {
                /* order does not divide p+1. */
                *is_supersingular = false;
                return true;
            }

            uint_mul3_64(order, order, primes[lower]);

            if (uint_sub3(&tmp, &four_sqrt_p, order)) { /* returns borrow */
                /* order > 4 sqrt(p), hence definitely supersingular */
                *is_supersingular = true;
                return true;
            }
        }

        /* inconclusive */
        return false;
    }

    size_t mid = lower + (upper - lower + 1) / 2;

    uint cl = uint_1, cu = uint_1;
    for (size_t i = lower; i < mid; ++i)
        uint_mul3_64(&cu, &cu, primes[i]);
    for (size_t i = mid; i < upper; ++i)
        uint_mul3_64(&cl, &cl, primes[i]);

    proj Q;

    xMUL(&Q, A, P, &cu,count);
    xMUL(P, A, P, &cl,count);

    /* start with the right half; bigger primes help more */
    return validate_rec(&Q, A, mid, upper, order, is_supersingular,count)
        || validate_rec(P, A, lower, mid, order, is_supersingular,count);
}

/* never accepts invalid keys. */
bool validate(public_key const *in,uint64_t* count)
{
    /* make sure the curve is nonsingular: A^2-4 != 0 */
    {
        uint dummy;
        if (!uint_sub3(&dummy, (uint *) &in->A, &p)) /* returns borrow */
            /* A >= p */
            return false;

        fp fp_pm2;
        fp_set(&fp_pm2, 2);
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = 2 */
            return false;

        fp_sub3(&fp_pm2, &fp_0, &fp_pm2);
        count[2]++;
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = -2 */
            return false;
    }

    const proj A = {in->A, fp_1};

    do {
        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        /* maximal 2-power in p+1 */
        xDBL(&P, &A, &P,count);
        xDBL(&P, &A, &P,count);

        bool is_supersingular;
        uint order = uint_1;

        if (validate_rec(&P, &A, 0, NUM_PRIMES, &order, &is_supersingular,count))
            return is_supersingular;

    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}



static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

/* (obviously) not constant time in the exponent */
static void fp_pow(fp *x, uint const *e,uint64_t* count)
{
    fp y = *x;
    *x = fp_1;//y=x,x=1
    for (size_t k = 0; k < LIMBS; ++k) {
        uint64_t t = e->c[k];
        for (size_t i = 0; i < 64; ++i, t >>= 1) {
            if (t & 1)
                fp_mul2(x, &y);
                count[0]++;
            fp_sq1(&y);
            count[1]++;
        }
    }
    
}

//if x is square ,compute y^2=x and return 1,else return 0
//has been checked right
void fp_square(fp *x, fp *y,uint64_t* count)
{
    fp tmp=*x;
    fp_pow(&tmp,&p_add_1_div_4,count);
    memcpy(y,&tmp,sizeof(fp));
}

//恢复y值
//checked right
void fp_compute_y(proj *Q,fp *A,fp *y,uint64_t* count)
{
    fp x,z;
    z=Q->z;
    fp_inv(&z);
    count[3]=count[3]+1;
    fp_mul3(&x,&z,&Q->x);
    count[0]=count[0]+1;
    fp rhs;
    montgomery_rhs(&rhs,A,&x,count);
    fp_square(&rhs,y,count);
}


//checked right
void montgomery_rhs(fp *rhs, fp const *A, fp const *x,uint64_t* count)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);

    count[0]=count[0]+2;
    count[1]=count[1]+1;
    count[2]=count[2]+2;
}

//checked right
// void finding_point(fp *A,proj *P,int i,uint64_t *count)
// {
//     fp rhs;
//     proj R;
// }
void finding_point(fp *A, proj *P, int i,uint64_t* count)
{
    fp rhs;
    proj R;
    int k = 1;
    uint cof;
    uint_set(&cof,4);
    proj A1;
    A1.x=*A;
    A1.z=fp_1;
    for(int m=0;m<NUM_PRIMES;m++)
    {
        if(m!=i) uint_mul3_64(&cof,&cof,primes[m]);
        else uint_mul3_64(&cof,&cof,1);
    }//cof=p-1/l_i
    do
    {
        do{
            fp_set(&R.x,k);
            R.z=fp_1;
            montgomery_rhs(&rhs,A,&R.x,count);
            if(fp_issquare(&rhs)) break;//确保是Fp上的点
            k++;
        }while(1);
        xMUL(P,&A1,&R,&cof,count);
        if((memcmp(&P->z,&fp_0,sizeof(fp))!=0)) break;//确保该点的阶为l_i
        k++;
    } while (1);
}

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

void point_add(proj *P,proj *Q,proj *R,fp *A,uint64_t* count)
{
    fp Y1,Y2;
    fp_compute_y(P,A,&Y1,count);
    fp_compute_y(Q,A,&Y2,count);

    fp a,b,c,d,e,f,g;
    a=P->z;
    b=Q->z;
    fp_inv(&a);
    fp_inv(&b);
    fp_mul3(&c,&P->x,&a);
    fp_mul3(&d,&Q->x,&b);

    fp_sub3(&e,&Y2,&Y1);//e=Y2-Y1
    fp_sub3(&f,&d,&c);//f=X2-X1
    fp_add3(&g,&c,&d);//g=X2+X1
    fp_inv(&f);//f=(X2-X1)^(-1)
    fp_mul2(&e,&f);//e=(Y2-Y1)*(X2-X1)^(-1)
    fp_sq1(&e);//e=[(Y2-Y1)*(X2-X1)^(-1)]^2
    fp_sub3(&R->x,&e,A);
    fp_sub2(&R->x,&g);
    R->z=fp_1;

    count[0]=count[0]+3;
    count[1]=count[1]+1;
    count[2]=count[2]+5;
    count[4]=count[4]+1;
}

/* totally not constant-time. */
void action(public_key *out, public_key const *in, private_key const *priv,uint64_t* count)
{
    uint k[2];
    uint_set(&k[0], 4); /* maximal 2-power in p+1 */
    uint_set(&k[1], 4); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES];

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int8_t t = (int8_t) priv->e[i] ;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    int Count = 0;

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x,count);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;
        
        Count ++;

        //uint64_t T = rdtsc();

        xMUL(&P, &A, &P, &k[sign],count);

        /*printf("%d , %d \n",count , rdtsc() - T);
        for(int l=0; l<NUM_PRIMES ; l++){
            printf("%2d", e[0][l]-e[1][l]);
        }
        printf("\n");*/

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof,count);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {                    

                    // T = rdtsc();

                    xISOG(&A, &P, &K, primes[i],0,count);

                    //printf("   i:%2d cyc:%d cyc/p:%d \n", i, rdtsc() - T, (rdtsc() - T)/primes[i]);

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        count[0]++;
        count[3]++;
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    //printf("\n");

    out->A = A.x;

}

static uint_print(uint* c)
{
    for(int i=7;i>=0;i--)
    {
        printf("%08X",(uint32_t)(c->c[i]>>32));
        printf("%08X",(uint32_t)(c->c[i]));
    }
    printf("\n\n");
}

void  action_one(public_key *out,public_key* in,proj* S,uint64_t* count)
{
	int check=0;
	proj A={in->A,fp_1};
    int i,j;
    proj R;
    uint a;
    {
        proj K;
        finding_point(&A.x,&K,42,count);
        xISOG(&A,&K,&K,primes[42],check,count);//计算第l_i次同源，并将R1映射到新曲线上
        fp_inv(&A.z);
        fp_mul2(&A.x,&A.z);
        count[0]=count[0]+1;
        count[3]=count[3]+1;
        A.z=fp_1;
        finding_point(&A.x,&R,42,count);
    }
    for(i=43;i<74;i++)
    {
        proj K;
        finding_point(&A.x,&K,i,count);
        xISOG(&A,&R,&K,primes[i],check,count);
        fp_inv(&A.z);
        fp_mul2(&A.x,&A.z);
        count[0]=count[0]+1;
        count[3]=count[3]+1;
        A.z=fp_1;
        proj Adding;
        finding_point(&A.x,&Adding,i,count);
        proj tmp;
        point_add(&R,&Adding,&tmp,&A.x,count);
        R.x=tmp.x;
        R.z=tmp.z;
    }
    S->x=R.x;
    S->z=R.z;
    out->A=A.x;

    //check S is of order wanted
    // uint cof;
    // uint_set(&cof,primes[42]);
    // for(int i=43;i<NUM_PRIMES;i++)
    // {
    //     uint_mul3_64(&cof,&cof,primes[i]);
    // }
    // proj Check;
    // xMUL(&Check,&A.x,S,&cof,count);
    // fp_dec(&a,&Check.z);
    // uint_print(&a);
}




/* includes public-key validation. */
bool csidh(public_key *out, public_key const *in, private_key const *priv,uint64_t* count)
{
    if (!validate(in,count)) {
        fp_random(&out->A);
        return false;
    }
    action(out, in, priv,count);
    return true;
}

