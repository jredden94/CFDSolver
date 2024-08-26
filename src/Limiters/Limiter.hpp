#pragma once

#include <functional>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>

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
        const double GetEpsilon(void) const;
        const Type LimiterType(void) const;

        function<double(const double& proj_grad, const double& delta, const double& eps)> Limit;

    private:
        Type type;
        double static VanAlbada(const double &proj_grad, const double &delta, const double &eps);
        double static Venkat(const double &proj_grad, const double &delta, const double &eps);
        double static NR3(const double &proj_grad, const double &delta, const double &eps);
        double static NR4(const double &proj_grad, const double &delta, const double &eps);
        double static NR5(const double &proj_grad, const double &delta, const double &eps);

        double venkat_eps = 0.5;
        static constexpr double machine_eps = numeric_limits<double>::epsilon();
        static constexpr double ref_len = 1.0;
        double eps;
};

