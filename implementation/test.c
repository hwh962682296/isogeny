#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "csifish.h"
#include "stdint.h"
#include <time.h>

#define KEYS 1
#define SIGNATURES_PER_KEY 100

static inline
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#define TIC printf("\n"); uint64_t cl = rdtsc();
#define TOC(A) printf("%s cycles = %lu \n",#A ,rdtsc() - cl); cl = rdtsc();


static uint_print(uint* c)
{
    for(int i=7;i>=0;i--)
    {
        printf("%08X",(uint32_t)(c->c[i]>>32));
        printf("%08X",(uint32_t)(c->c[i]));
    }
    printf("\n\n");
}

static void test_compute_y()
{
	uint a;
	uint64_t count[4]={0};
	proj P;
	P.z=fp_1;
	do{
		fp_random(&P.x);
		fp rhs;
		montgomery_rhs(&rhs,&base.A,&P.x,count);
		if(fp_issquare(&rhs)==1) break;
	}while(1);
	fp Y;
	fp_compute_y(&P,&base.A,&Y,count);
	fp Y_square;
	fp_sq2(&Y_square,&Y);
	fp X_rhs;
	fp_sq2(&X_rhs,&P.x);
	fp_mul2(&X_rhs,&P.x);
	fp_add2(&X_rhs,&P.x);
	int flag=memcmp(&X_rhs,&Y_square,sizeof(fp));
	if(flag==0) printf("Test compute Y pass\n");
	else printf("Test compute Y failed\n");
}

static void test_finding_point()
{
	uint64_t count[4]={0};
	proj A={base.A,fp_1};
	for(int i=0;i<NUM_PRIMES;i++)
	{
		proj P;
		finding_point(&base.A,&P,i,count);
		uint cof;
		uint_set(&cof,primes[i]);
		proj R;
		xMUL(&R,&A,&P,&cof,count);
		if(memcmp(&R.z,&fp_0,sizeof(fp)!=0))
		{
			printf("test_finding_point Faile\n");
			exit(0);
		}
	}
	printf("Test finding_point Pass\n");
}


static void test_point_add()
{
	uint64_t count[4]={0};
	proj P,Q1,Q2,R;
	uint a;
	P.z=fp_1;
	proj A={base.A,fp_1};
	uint cof1,cof2;
	uint_set(&cof1,5);
	uint_set(&cof2,4);
	fp_set(&P.x,4);
	xMUL(&Q1,&A,&P,&cof1,count);//Q1=5P
	xMUL(&Q2,&A,&P,&cof2,count);//Q2=4p
	fp_inv(&Q1.z);
	fp_mul2(&Q1.x,&Q1.z);
	fp_dec(&a,&Q1.x);
	uint_print(&a);
	point_add(&P,&Q2,&R,&base.A,count);//R=Q2+P
	fp_dec(&a,&R.x);
	uint_print(&a);
	if(memcmp(&Q1.x,&R.x,sizeof(fp))==0)
	{
		printf("Test point add PASS!\n");
		return 1;
	}
	else
	{
		printf("Test point add failed\n");
		return 0;
	}
}

static int test_order_add()
{
	uint64_t count[4]={0};
	proj P,Q,R;
	proj A={base.A,fp_1};
	finding_point(&base.A,&P,50,count);
	finding_point(&base.A,&Q,60,count);

	point_add(&P,&Q,&R,&base.A,count);

	uint cof;
	uint_set(&cof,primes[50]);
	{
		uint cof2;
		uint_set(&cof2,primes[50]);
		proj tmp;
		xMUL(&tmp,&A,&R,&cof2,count);
		uint a;
		fp_dec(&a,&tmp.z);
		uint_print(&a);
	}
	{
		uint cof3;
		uint_set(&cof3,primes[60]);
		proj tmp;
		xMUL(&tmp,&A,&R,&cof3,count);
		uint a;
		fp_dec(&a,&tmp.z);
		uint_print(&a);
	}

	uint_mul3_64(&cof,&cof,primes[60]);

	proj S;
	xMUL(&S,&A,&R,&cof,count);
	// uint a;
	// fp_dec(&a,&S.z);
	// uint_print(&a);
	if(memcmp(&S.z,&fp_0,sizeof(fp))==0)
	{
		printf("Order test right\n");
		return 1;
	}
	else
	{
		printf("Order test failed\n");
		return 0;
	}
}

static int check_point_eql(proj *P,proj *Q)
{
    fp a,b;
    fp_mul3(&a,&P->x,&Q->z);
    fp_mul3(&b,&P->z,&Q->x);

    return (memcmp(&a,&b,sizeof(fp)));
}



static void test_Fles()
{
	uint64_t count[4]={0};
	
}

static void test_CRT()
{
    int crs[32]={
        184,125,96,156,18,180,
        111,111,53,69,42,0,24,114,231,227,39,201,286,305,194,243,
        44,329,10,105,293,322,272,328,98,287
    };
    unsigned char message[32]={0};
    CRT(crs,message);
    for(int i=0;i<32;i++)
    {
        printf("0x%0x,",message[i]);
    }
}

static void test_performance()
{

	uint64_t count1[4]={0};
	uint64_t count2[4]={0};
	uint64_t count3[4]={0};
	uint32_t tag=0;

	unsigned char *pk=aligned_alloc(64,sizeof(uint));
	unsigned char *sk=aligned_alloc(64,16);
	for(int i=0;i<10000;i++)
	{
		Sims_keygen(sk,pk,count1);
	}
	unsigned char message[32]={0};
	for(int i=0;i<32;i++)
	{
		message[i]=rand();
		printf("%0X",message[i]);
	}
	unsigned char Message[32]={0};
	unsigned char *c=aligned_alloc(64,2*sizeof(uint));
	for(int i=0;i<10000;i++)
	{
		Sims_enc(pk,message,c,&tag,count2);
		printf("\n");
		Sims_dec(sk,c,Message,&tag,count3);
	}

	for(int i=0;i<4;i++)
	{
		count1[i]=count1[i]/10000;
		printf("%d\n",count1[i]);
	}
	printf("\n");

	for(int i=0;i<4;i++)
	{
		count2[i]=count2[i]/10000;
		printf("%d\n",count2[i]);
	}
	printf("\n");

	for(int i=0;i<4;i++)
	{
		count3[i]=count3[i]/10000;
		printf("%d\n",count3[i]);
	}
}

int main()
{


	uint64_t count1[4]={0};
	uint64_t count2[4]={0};
	uint64_t count3[4]={0};
	uint32_t tag=0;

	unsigned char *pk=aligned_alloc(64,sizeof(uint));
	unsigned char *sk=aligned_alloc(64,16);
	for(int i=0;i<1;i++)
	{
		Sims_keygen(sk,pk,count1);
	}
	Sims_keygen(sk,pk,count1);
	unsigned char message[32]={0};
	srand(time(0));
	for(int i=0;i<32;i++)
	{
		message[i]=rand();
		printf("%0X",message[i]);
	}
	unsigned char Message[32]={0};
	unsigned char *c=aligned_alloc(64,2*sizeof(uint));
	Sims_enc(pk,message,c,&tag,count2);
	printf("\n");
	Sims_dec(sk,c,Message,&tag,count3);

	int flag=1;
	for(int i=0;i<32;i++)
	{
		if(message[i]!=Message[i]) flag = 0;
		else flag = 1;
	}
	if(flag==1) printf("FleS is right!\n");
	else printf("Failed\n");

}

// static void message_2_uint(unsigned char* m,uint* k)
// {
// 	uint64_t* tmp=m;
	
// 	for(int i=0;i<4;i++)
// 	{
// 		k->c[i]=tmp[i];
// 	}
// 	k->c[4]=0;
// 	k->c[5]=0;
// 	k->c[6]=0;
// 	k->c[7]=0;
// }

// int main()
// {
// 	unsigned char m[32]={
// 		0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xff,
// 		0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xff,
// 		0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xff,
// 		0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xff,
// 	};
// 	uint c;
// 	message_2_uint(m,&c);
// 		for(int i=7;i>=0;i--)
//     {
//         printf("%08X",(uint32_t)(c.c[i]>>32));
//         printf("%08X",(uint32_t)(c.c[i]));
//     }
// }
