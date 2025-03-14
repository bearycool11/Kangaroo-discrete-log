#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <stdint.h>

static uint64_t modexp_cpp(uint64_t base, uint64_t exp, uint64_t mod)
{
    uint64_t result = 1ULL;
    uint64_t cur = base % mod;
    uint64_t e   = exp;
    while (e > 0)
    {
        if (e & 1ULL)
        {
            __uint128_t tmp = ( __uint128_t)result * cur;
            result = (uint64_t)(tmp % mod);
        }
        __uint128_t tmp2 = ( __uint128_t)cur * cur;
        cur = (uint64_t)(tmp2 % mod);
        e >>= 1ULL;
    }
    return result;
}

class PollardKangaroo
{
public:
    PollardKangaroo()
    {
        partition_count_ = 32;
        jump_table_ = nullptr;
        exp_table_ = nullptr;
    }

    ~PollardKangaroo()
    {
        if (jump_table_)
            delete[] jump_table_;
        if (exp_table_)
            delete[] exp_table_;
    }

    long long solve(uint64_t g, uint64_t h, uint64_t p, long long a, long long b)
    {
        if (a > b) return -1;

        long long interval_size = b - a;
        if (interval_size < 0) return -1;

        // power_gb = g^b mod p
        uint64_t power_gb = modexp_cpp(g, (uint64_t)b, p);

        // Tame start
        uint64_t T = power_gb;
        long long T_exp = b;

        // Wild start
        uint64_t W = h % p;
        long long W_exp = 0LL;

        long long max_jump = (interval_size / (2LL * partition_count_)) + 1LL;

        // free old if exist
        if (jump_table_) delete[] jump_table_;
        if (exp_table_) delete[] exp_table_;
        jump_table_ = new uint64_t[partition_count_];
        exp_table_  = new long long[partition_count_];

        std::srand((unsigned)std::time(nullptr));

        for (int i=0; i< partition_count_; i++)
        {
            long long s_i = 1 + (std::rand() % (int)max_jump);
            jump_table_[i] = modexp_cpp(g, (uint64_t)s_i, p);
            exp_table_[i]  = s_i;
        }

        auto partition = [this](uint64_t x) {
            return (int)(x & (this->partition_count_ - 1));
        };

        long long loop_bound = (long long)(2.0 * std::sqrt((long double)interval_size)) + 1000LL;
        if (loop_bound < 0) loop_bound = 1000LL;

        for (long long iteration=0; iteration < loop_bound; iteration++)
        {
            // Tame
            int idx_t = partition(T);
            __uint128_t tmpT = ( __uint128_t)T * jump_table_[idx_t];
            T = (uint64_t)(tmpT % p);
            T_exp -= exp_table_[idx_t];

            // Wild
            int idx_w = partition(W);
            __uint128_t tmpW = ( __uint128_t)W * jump_table_[idx_w];
            W = (uint64_t)(tmpW % p);
            W_exp += exp_table_[idx_w];

            if (T == W)
            {
                // collision
                long long candidate_x = (long long)(b - T_exp);
                if (candidate_x >= a && candidate_x <= b)
                {
                    uint64_t check = modexp_cpp(g, (uint64_t)candidate_x, p);
                    if (check == (h % p))
                    {
                        return candidate_x;
                    }
                }
                return -1;
            }
        }
        return -1; // not found
    }

private:
    int partition_count_;
    uint64_t* jump_table_;
    long long* exp_table_;
};

int main()
{
    // example usage
    uint64_t p = 467ULL;
    uint64_t g = 2ULL;
    long long a = 20LL;
    long long b = 100LL;

    long long x_true = 63LL;
    // compute h = g^x mod p
    uint64_t h = 1ULL;
    for (int i=0; i < x_true; i++)
        h = (h * g) % p;

    std::cout << "Testing x_true=" << x_true << " => h=" << h << std::endl;

    PollardKangaroo kangaroo;
    long long found_x = kangaroo.solve(g, h, p, a, b);
    if (found_x >= 0)
        std::cout << "Discrete log found: x = " << found_x << std::endl;
    else
        std::cout << "No discrete log found or out of range." << std::endl;

    return 0;
}
