#ifndef CSIFISH_H
#define CSIFISH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csidh.h"
#include "merkletree.h"
#include "stdint.h"
#include "parameters.h"
#include "classgroup.h"
#include "fp.h"
#include "uint.h"
#include "mont.h"

#ifdef MERKLEIZE_PK

	#define PK_BYTES SEED_BYTES+SEED_BYTES

	#define SK_SEED(sk) (sk)
	#define SK_MERKLE_KEY(sk) (SK_SEED(sk) + SEED_BYTES)
	#define SK_TREE(sk) (SK_MERKLE_KEY(sk) + SEED_BYTES)
	#define SK_BYTES (2*SEED_BYTES + ((2*PKS-1)*SEED_BYTES))

	#define SIG_CURVES(sig) (sig) 
	#define SIG_RESPONSES(sig) (SIG_CURVES(sig) + sizeof(uint[ROUNDS]))
	#define SIG_TREE_FILLING(sig) (SIG_RESPONSES(sig) + 33*ROUNDS)
	#define SIG_BYTES (SIG_TREE_FILLING(0) + 10000)

#else
	
	#define PK_CURVES(pk) (pk)
	#define PK_BYTES sizeof(uint[PKS])

	#define SK_SEED(sk) (sk)
	#define SK_BYTES SEED_BYTES 

	#define SIG_HASH(sig) (sig)
	#define SIG_RESPONSES(sig) (SIG_HASH(sig) + HASH_BYTES)
	#define SIG_BYTES (SIG_RESPONSES(0) + 33*ROUNDS)

#endif

void FleS_keygen(unsigned char* sk,unsigned char* pk,uint64_t* count);
void FleS_enc(unsigned char* pk,unsigned char* m,unsigned char* c,uint32_t* tag,uint64_t* count);
void FleS_dec(unsigned char* sk,unsigned char* c,unsigned char* m,uint32_t* tag,uint64_t* count);

#endif