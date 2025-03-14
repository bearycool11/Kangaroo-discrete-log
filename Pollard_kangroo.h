#ifndef POLLARD_KANGAROO_H
#define POLLARD_KANGAROO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Pollard's Kangaroo Algorithm for discrete log g^x = h mod p, with x in [a, b].
 * 
 * \param g           Base/generator (unsigned long long)
 * \param h           Target h = g^x mod p (unsigned long long)
 * \param p           Prime modulus (unsigned long long)
 * \param a           Lower bound for x
 * \param b           Upper bound for x
 * \return            x if found, or -1 if not found or out of range
 */
long long pollard_kangaroo(uint64_t g, uint64_t h, uint64_t p,
                           long long a, long long b);

#ifdef __cplusplus
}
#endif

#endif // POLLARD_KANGAROO_H
