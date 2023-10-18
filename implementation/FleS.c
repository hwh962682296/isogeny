#include "csifish.h"

void print_uint(uint x){
	for(int i=0 ; i<LIMBS; i++){
		printf("%lu ", x.c[i] );
	}
	printf("\n");
}

const private_key priv_one={{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};

void Sims_keygen(unsigned char* sk,unsigned char* pk,uint64_t* count)
{
	RAND_bytes(sk,SEED_BYTES);

	init_classgroup();

	uint* CURVE = PK_CURVES(pk);

	private_key priv;
	public_key out;
	sample_from_classgroup_with_seed(sk,priv.e);//sk做种子生成同源
	
	action(&out,&base,&priv,count);
	fp_dec(CURVE,&out.A);

	clear_classgroup();
}

#ifndef f_E
#define f_E(E,x,x1)			\
{							\
	x1[0]=x[0]^E[0];		\
	x1[1]=x[1]^E[1];		\
	x1[2]=x[2]^E[2];		\
	x1[3]=x[3]^E[3];		\
	x1[4]=x[4]^E[4];		\
	x1[5]=x[5]^E[5];		\
	x1[6]=x[6]^E[6];		\
	x1[7]=x[7]^E[7];		\
}
#endif

#ifndef g_E
#define g_E(E,x,x1)			\
{							\
	x[0]=x1[0]^E[0];		\
	x[1]=x1[1]^E[1];		\
	x[2]=x1[2]^E[2];		\
	x[3]=x1[3]^E[3];		\
	x[4]=x1[4]^E[4];		\
	x[5]=x1[5]^E[5];		\
	x[6]=x1[6]^E[6];		\
	x[7]=x1[7]^E[7];		\
}
#endif

static void message_2_uint(unsigned char* m,uint* k)
{
	unsigned char TMP[32]={0};
	for(int i=0;i<32;i++)
	{
		TMP[i]=m[31-i];
	}
	uint64_t* tmp=TMP;
	
	for(int i=0;i<4;i++)
	{
		k->c[i]=tmp[i];
	}
	k->c[4]=0;
	k->c[5]=0;
	k->c[6]=0;
	k->c[7]=0;
}

static void compute_tag(uint32_t* tag,unsigned char* message)
{

	mpz_t a;
	mpz_t p;
	mpz_t p_2;
	mpz_t r;
	mpz_init(a);
	mpz_init(p);
	mpz_init(r);

	mpz_import(a,32,1,1,0,0,message);
	printf("\n");
	gmp_printf("%Zd\n",a);
	for(int i=0;i<32;i++)
	{
		mpz_set_ui(p,primes[i+42]);
		mpz_set_ui(p_2,(primes[i+42]-1)/2);
		mpz_mod(r,a,p);
		if(mpz_cmp(r,p_2)==1) *tag = *tag+(1<<i);
		else *tag = *tag+0;
	}

	mpz_clear(a);
	mpz_clear(p);
	mpz_clear(r);

}

void Sims_enc(unsigned char* pk,unsigned char* m, unsigned char* c,uint32_t* tag,uint64_t* count)
{

	compute_tag(tag,m);
	unsigned char* seed = malloc(SEED_BYTES);
	RAND_bytes(seed,SEED_BYTES);

	init_classgroup();
	private_key priv,priv1;
	sample_from_classgroup_with_seed(seed,priv.e);
	for(int i=0;i<42;i++)
	{
		priv1.e[i]=priv.e[i];
	}
	for(int i=42;i<NUM_PRIMES;i++)
	{
		priv1.e[i]=priv.e[i]-1;
	}
	public_key start,tmp,end,E_b;
	proj P,Q;
	fp_enc(&(start.A),(uint*)pk);//start为Alice公钥E_a

	action(&tmp,&start,&priv1,count);//tmp为E_mid
	action(&E_b,&base,&priv,count);
	uint E_b_tmp;
	fp_dec(&E_b_tmp,&E_b.A);
	action_one(&end,&tmp,&P,count);
	proj A={end.A,fp_1};


	uint scale;
	message_2_uint(m,&scale);
    // for(int i=7;i>=0;i--)
    // {
    //     printf("%08X",(uint32_t)(scale.c[i]>>32));
    //     printf("%08X",(uint32_t)(scale.c[i]));
    // }
    // printf("\n\n");

	xMUL(&Q,&A,&P,&scale,count);//消息m嵌入密文

	fp Q_c;
	fp_inv(&Q.z);
	fp_mul2(&Q.x,&Q.z);

	count[0]=count[0]+1;
	count[3]=count[3]+1;
	f_E(end.A.c,Q.x.c,Q_c.c);

	memcpy(c,(unsigned char*)&E_b.A,sizeof(fp));
	memcpy(c+sizeof(fp),(unsigned char*)&Q_c,sizeof(fp));

	clear_classgroup();
}


void Sims_dec(unsigned char* sk,unsigned char* c,unsigned char* m,uint32_t* tag,uint64_t* count)
{
	
	init_classgroup();
	private_key priv,priv1;
	sample_from_classgroup_with_seed(sk,priv.e);
	for(int i=0;i<42;i++)
	{
		priv1.e[i]=priv.e[i];
	}
	for(int i=42;i<NUM_PRIMES;i++)
	{
		priv1.e[i]=priv.e[i]-1;
	}
	
	public_key start,tmp,end;
	memcpy(&(start.A),c,sizeof(fp));


	action(&tmp,&start,&priv1,count);

	proj P;
	action_one(&end,&tmp,&P,count);//计算E_ab同时计算点P_E_4

	proj A={end.A,fp_1};
	
	proj Q;
	Q.z=fp_1;
	fp Q_x_c;
	memcpy(&Q_x_c,c+sizeof(fp),sizeof(fp));
	g_E(end.A.c,Q.x.c,Q_x_c.c);

	Pohlig_Hellman(&Q,&P,&A,m,tag,count);//利用Pohlig-Hellman算法解出消息m

}
