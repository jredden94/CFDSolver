#include "Limiter.hpp"
#include "../Common/Config.hpp"

Limiter::Limiter() { 
    SetLimiter(Type::VanAlbada);
    Config *config = &Config::GetConfig();
    venkat_eps = config->GetLimiterCoefficient();
}
Limiter::Limiter(Type type) : type(type) { 
    SetLimiter(type);
}
Limiter::~Limiter() { }

const zdouble Limiter::GetEpsilon(void) const { return eps; }

void Limiter::SetLimiter(const Type &type) {
    this->type = type;
    eps = machine_eps;
    if (type == Type::VanAlbada) Limit = VanAlbada;
    else if (type == Type::Venkatakrishnan) {
        zdouble eps1 = fabs(ref_len * venkat_eps);
        eps = max(pow(eps1,3), machine_eps);
        Limit = Venkat;
    }
    else if (type == Type::Nishikawa_R3) {
        zdouble eps1 = fabs(ref_len * venkat_eps);
        eps = max(pow(eps1,4), machine_eps);
        Limit = NR3;
    }
    else if (type == Type::Nishikawa_R4) {
        zdouble eps1 = fabs(ref_len * venkat_eps);
        eps = max(pow(eps1,5), machine_eps);
        Limit = NR4;
    }
    else if (type == Type::Nishikawa_R5) {
        zdouble eps1 = fabs(ref_len * venkat_eps);
        eps = max(pow(eps1,6), machine_eps);
        Limit = NR5;
    }
}

const Limiter::Type Limiter::LimiterType() const { return type; }

zdouble Limiter::VanAlbada(const zdouble& proj_grad, const zdouble &delta, const zdouble &eps) {
    return delta * (2.0 * proj_grad + delta) / (4.0 * proj_grad * proj_grad + delta * delta + eps);
}

zdouble Limiter::Venkat(const zdouble& proj_grad, const zdouble &delta, const zdouble &eps) {
    zdouble y = delta * (delta + proj_grad) + eps;
    return (y + delta * proj_grad) / (y + 2.0 * proj_grad * proj_grad);
}

zdouble Limiter::NR3(const zdouble& proj_grad, const zdouble &delta, const zdouble &eps) {
    zdouble dp = fabs(delta);
    zdouble dm = fabs(proj_grad);
    if (dp > (2.0 * dm)) return 1.0;
    zdouble y = dp * dp * dp + eps;
    zdouble s3 = 4.0 * dm * dm;
    return (y + dp * s3) / (y + dm * (delta * delta + s3));
}

zdouble Limiter::NR4(const zdouble& proj_grad, const zdouble &delta, const zdouble &eps) {
    zdouble dp = fabs(delta);
    zdouble dm = fabs(proj_grad);
    if (dp > (2.0 * dm)) return 1.0;
    zdouble y = pow(dp, 4.0) + eps;
    zdouble s4 = 2.0 * dm * (dp * dp - 2.0 * dm * (dp - 2.0 * dm));
    return (y + dp * s4) / (y + dm * (pow(delta,3.0) + s4));
}

zdouble Limiter::NR5(const zdouble& proj_grad, const zdouble &delta, const zdouble &eps) {
    zdouble dp = fabs(delta);
    zdouble dm = fabs(proj_grad);
    if (dp > (2.0 * dm)) return 1.0;
    zdouble y = pow(dp, 5.0) + eps;
    zdouble s5 = 8.0 * dm * dm * (dp * dp - 2.0 * dm * (dp - dm));
    return (y + dp * s5) / (y + dm * (pow(delta,4.0) + s5));
}
