#include "pollard_kangaroo.h"
#include <stdlib.h>   // for rand, srand
#include <time.h>     // for time
#include <math.h>     // for sqrt
#include <stdio.h>    // for printf

/* 
 * We'll define a helper function for pow(base, exp) % mod, 
 * but note that C99 doesn't guarantee 64-bit safe math 
 * for big mod. For real usage with large p, you'd want a 
 * big-integer library. For demonstration, we'll assume p < 2^63.
 */
static uint64_t modexp(uint64_t base, uint64_t exp, uint64_t mod)
{
    uint64_t result = 1ULL;
    uint64_t cur = base % mod;
    uint64_t e   = exp;

    while (e > 0) {
        if (e & 1ULL) {
            __uint128_t tmp = ( __uint128_t)result * cur;
            result = (uint64_t)(tmp % mod);
        }
        __uint128_t tmp2 = ( __uint128_t)cur * cur;
        cur = (uint64_t)(tmp2 % mod);
        e >>= 1ULL;
    }

    return result;
}

static int partition_count = 32;  // controls how many partition sets
static uint64_t *jump_table = NULL;
static long long *exp_table = NULL;

// simple partition function using lower bits
static int partition(uint64_t x)
{
    // Take last 5 bits for range [0..31]
    return (int)(x & (partition_count - 1));
}

long long pollard_kangaroo(uint64_t g, uint64_t h, uint64_t p,
                           long long a, long long b)
{
    if (a > b) return -1; // invalid range

    // interval_size = b - a
    long long interval_size = b - a;
    if (interval_size < 0) return -1;

    // g^b mod p
    uint64_t power_gb = modexp(g, (uint64_t)b, p);

    // Tame kangaroo starts at T = g^b, exponent track = b
    uint64_t T = power_gb;
    long long T_exp = b;

    // Wild kangaroo starts at W = h, exponent track = 0
    uint64_t W = h % p;
    long long W_exp = 0LL;

    // Generate jump table for random steps in [1..someMax]
    // We'll do a dynamic approach based on interval_size
    // approximate max jump
    long long max_jump = (interval_size / (2LL * partition_count)) + 1LL;

    // if already allocated, free
    if (jump_table) {
        free(jump_table);
        jump_table = NULL;
    }
    if (exp_table) {
        free(exp_table);
        exp_table = NULL;
    }

    jump_table = (uint64_t*)malloc(sizeof(uint64_t) * partition_count);
    exp_table  = (long long*)malloc(sizeof(long long) * partition_count);

    // seed random
    srand((unsigned)time(NULL));

    for (int i = 0; i < partition_count; i++) {
        // random [1..max_jump]
        long long s_i = 1 + (rand() % (int)(max_jump));
        jump_table[i] = modexp(g, (uint64_t)s_i, p);
        exp_table[i]  = s_i;
    }

    // loop bound ~ 2*sqrt(interval_size) + 1000
    long long loop_bound = (long long)(2.0*sqrt((long double)interval_size)) + 1000LL;
    if (loop_bound < 0) loop_bound = 1000LL; // fallback

    for (long long iteration = 0; iteration < loop_bound; iteration++) {
        // Tame step
        int idx_t = partition(T);
        T = ( ( __uint128_t)T * jump_table[idx_t]) % p;
        T_exp -= exp_table[idx_t];

        // Wild step
        int idx_w = partition(W);
        W = ( ( __uint128_t)W * jump_table[idx_w]) % p;
        W_exp += exp_table[idx_w];

        if (T == W) {
            // collision
            long long candidate_x = (long long)(b - T_exp);
            if (candidate_x >= a && candidate_x <= b) {
                // verify
                uint64_t check = modexp(g, (uint64_t)candidate_x, p);
                if (check == (h % p)) {
                    // success
                    free(jump_table);
                    free(exp_table);
                    jump_table = NULL;
                    exp_table = NULL;
                    return candidate_x;
                }
            }
            // if it didn't verify or out-of-range, no solution here
            free(jump_table);
            free(exp_table);
            jump_table = NULL;
            exp_table = NULL;
            return -1;
        }
    }

    // not found
    free(jump_table);
    free(exp_table);
    jump_table = NULL;
    exp_table = NULL;
    return -1;
}
