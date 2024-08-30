#pragma once

#include <functional>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>
#include "../Common/AD.hpp"

class Config;

using namespace std;

class Limiter {
    public:
        enum class Type { None, VanAlbada, Venkatakrishnan, Nishikawa_R3, Nishikawa_R4, Nishikawa_R5 };

        Limiter();
        Limiter(const Type type);
        ~Limiter();
        Limiter(const Limiter&) = delete;
        Limiter operator=(const Limiter&) = delete;

        void SetLimiter(const Type &type);
        const zdouble GetEpsilon(void) const;
        const Type LimiterType(void) const;

        function<zdouble(const zdouble& proj_grad, const zdouble& delta, const zdouble& eps)> Limit;

    private:
        Type type;
        zdouble static VanAlbada(const zdouble &proj_grad, const zdouble &delta, const zdouble &eps);
        zdouble static Venkat(const zdouble &proj_grad, const zdouble &delta, const zdouble &eps);
        zdouble static NR3(const zdouble &proj_grad, const zdouble &delta, const zdouble &eps);
        zdouble static NR4(const zdouble &proj_grad, const zdouble &delta, const zdouble &eps);
        zdouble static NR5(const zdouble &proj_grad, const zdouble &delta, const zdouble &eps);

        zdouble venkat_eps = 0.5;
        zdouble machine_eps = numeric_limits<zdouble>::epsilon();
        zdouble ref_len = 1.0;
        zdouble eps;
};

