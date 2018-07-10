/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for FourQ's field arithmetic 
************************************************************************************/    

#include "../FourQ_api.h"
#include "../FourQ_internal.h"
#include "../FourQ_params.h"
#include "test_extras.h"
#include <stdio.h>
#include <string.h>

#include <time.h>
#include "blake2.h"


// Benchmark and test parameters 
#define BENCH_LOOPS       10000      // Number of iterations per bench
#define SHORT_BENCH_LOOPS 1000       // Number of iterations per bench (for expensive operations)
#define TEST_LOOPS        1000       // Number of iterations per test

int K = 100;


bool fp2_run()
{
    bool OK = true;
    int n, i;
    unsigned long long cycles, cycles1, cycles2;
    double time;
    time = 0.0;
    clock_t start, start2;
    clock_t end, end2;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking MRETA with K = %d : \n\n", K);

    point_t P;
    point_t RJ;
    unsigned char y[32] = {0x54, 0xa2, 0xf8, 0x03, 0x1d, 0x18, 0xac, 0x77, 0xd2, 0x53, 0x92, 0xf2, 0x80, 0xb4, 0xb1, 0x2f, 0xac, 0xf1, 0x29, 0x3f, 0x3a, 0xe6, 0x77, 0x7d, 0x74, 0x15, 0x67, 0x91, 0x99, 0x53, 0x69, 0xc5}; 
    digit_t* ydigit = (digit_t*)y;
    unsigned char Y[32], rj[32], zj[32], Rj[32], gamma[32*K], beta[32*K];
    digit_t* rjdigit = (digit_t*)rj;
//    unsigned char tobeHashed[1] = {0};
    unsigned char tobeHashed[2] = {0};
    unsigned char hashed[64] = {0};

    cycles = 0;
    for (n = 0; n<SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
 //       start = clock();

        ecc_mul_fixed((digit_t*)y, P);
        encode(P, Y);

        for(int i = 0; i<K; i++){
//            blake2b(hashed, tobeHashed, y, 64, 1, 32);
//            memcpy(rj, hashed, 32);
//            memcpy(zj, hashed+32, 32);
            tobeHashed[1] = 0;
            blake2s(rj, tobeHashed, y, 32, 2, 32);
            tobeHashed[1] = 1;
            blake2s(zj, tobeHashed, y, 32, 2, 32);

            ecc_mul_fixed((digit_t*)rj, P);
            encode(P, Rj);

//            blake2b(beta+32*i, Rj, NULL, 32, 32, 0);

            blake2s(beta+32*i, Rj, NULL, 32, 32, 0);

            for(int j = 0; j<32; j++){
                gamma[j + 32*i] = beta[j + 32*i]^zj[j];
            }

            tobeHashed[0] += 1;
        }
        
 //       end = clock();
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
        time = time + (double)(end-start);

        tobeHashed[0] = 0;
    }
    printf("  MRETA Key generation runs in .............. %8lld ", (cycles/SHORT_BENCH_LOOPS)); print_unit;
    printf("\n");
//    printf("%fus per key generation\n", ((double) (time * 1000)) / CLOCKS_PER_SEC / SHORT_BENCH_LOOPS * 1000);

    unsigned char msg[32] = {"a"};
    unsigned char msgCount[2] = {0};
    unsigned char cj[32], ej[32], sj[32], tobeHashed2[64];
    digit_t* sjdigit = (digit_t*)sj;
    digit_t* ejdigit = (digit_t*)ej;

    time = 0.0;
    cycles = 0;
    msgCount[0] = 1;
    for (n = 0; n<SHORT_BENCH_LOOPS*1000; n++)
    {
        cycles1 = cpucycles();
//        start = clock();
//        blake2b(hashed, msgCount, y, 64, 1, 32);
//        memcpy(rj, hashed, 32);
//        memcpy(zj, hashed+32, 32);

        msgCount[1] = 0;
        blake2s(rj, msgCount, y, 32, 2, 32);
        msgCount[1] = 1;
        blake2s(zj, msgCount, y, 32, 2, 32);

        for(int j = 0; j<32; j++){
            cj[j] = zj[j]^msg[j];
        }   

//        memcpy(tobeHashed2, cj, 32);
//        memcpy(tobeHashed2+32,msg,32);

//        blake2b(ej, cj, NULL, 32, 32, 0);

        blake2s(ej, cj, NULL, 32, 32, 0);

        

        to_Montgomery(ejdigit, ejdigit);
        to_Montgomery(ydigit, sjdigit);
        Montgomery_multiply_mod_order(sjdigit, ejdigit, sjdigit);
        from_Montgomery(sjdigit, sjdigit);
        subtract_mod_order(rjdigit, sjdigit, sjdigit);

//        end = clock();
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
        time = time + (double)(end-start);
    }
    printf("  MRETA Signature generation runs in ..... %8lld ", cycles/(SHORT_BENCH_LOOPS*1000)); print_unit;
    printf("\n");
//    printf("%fus per Signature generation\n", ((double) (time * 1000)) / CLOCKS_PER_SEC / SHORT_BENCH_LOOPS * 1000);

    printf("sj: ");
    for (int i = 0; i < 32; i++) {
      printf("%x, ", sj[i]);
    }
    printf("\n");
    printf("cj: ");
    for (int i = 0; i < 32; i++) {
      printf("%x, ", cj[i]);
    }
    printf("\n");

    unsigned char check[32];
    time = 0.0;
    cycles = 0;
    for (n = 0; n<SHORT_BENCH_LOOPS*100; n++)
    {
        

        cycles1 = cpucycles();
//        start = clock();

//        memcpy(tobeHashed2, cj, 32);
//        memcpy(tobeHashed2+32,msg,32);
//        blake2b(ej, cj, NULL, 32, 32, 0);

        blake2s(ej, cj, NULL, 32, 32, 0);

        decode(Y, P);

        ecc_mul_double(sjdigit, P, ejdigit, RJ);
        encode(RJ, rj);

//        blake2b(check, rj, NULL, 32, 32, 0);
        blake2s(check, rj, NULL, 32, 32, 0);
        //Check if the signature is verified
        if (memcmp(check,beta+32,32) == 0){
            //printf("SIGNATURE IS VERIFIED\n");
            for(int j = 0; j<32; j++){
                //Recover the message
                msg[j] = gamma[j]^cj[j+32]^check[j];
                //printf("%c", msg[j]);
            }   
        }
        else
            printf("SIGNATURE IS NOT VERIFIED\n");

//        end = clock();
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
        time = time + (double)(end-start);
    }
    printf("  MRETA Signature verification runs in . %8lld ", cycles/(SHORT_BENCH_LOOPS*100)); print_unit;
    printf("\n");
//    printf("%fus per Signature verification\n", ((double) (time * 1000)) / CLOCKS_PER_SEC / SHORT_BENCH_LOOPS * 1000);
    
    return OK;
}


int main()
{
    bool OK = true;

    unsigned char tempor2[16];
    uint8_t hash2[2] = {0};
    blake2b(tempor2, hash2, NULL, 16, 2, 0);

    OK = OK && fp2_run();      // Benchmark quadratic extension field operations using p = 2^127-1
    
    return OK;
}
