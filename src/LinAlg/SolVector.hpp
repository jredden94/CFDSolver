#pragma once

#include <cassert>
#include <cmath>
#include "../Common/Config.hpp"
#include "../Common/AD.hpp"

class SolVector {
    public:
        SolVector();
        SolVector(const SolVector&);
        ~SolVector();

        void operator=(const SolVector&);
        zdouble* operator[](const unsigned long&);
        const zdouble* const operator[](const unsigned long&) const;
        void operator+=(const SolVector &vec);
        void operator-=(const SolVector &vec);
        void operator*=(const zdouble scalar);
        void operator/=(const zdouble scalar);
        void PlusEqualScalarMult(const SolVector &vec, const zdouble scalar);
        void MinusEqualScalarMult(const SolVector &vec, const zdouble scalar);

        void SetBlock(const unsigned long&, const vector<zdouble>&, const zdouble& = 1);
        void SetBlock(const unsigned long&, const zdouble*, const zdouble& = 1);
        void AddBlock(const unsigned long&, const vector<zdouble>&, const zdouble& = 1);
        void AddBlock(const unsigned long&, const zdouble*, const zdouble& = 1);
        void SubtractBlock(const unsigned long&, const vector<zdouble>&, const zdouble& = 1);
        void SubtractBlock(const unsigned long&, const zdouble*, const zdouble& = 1);
        const zdouble* const GetBlock(const unsigned long&) const;
        zdouble* GetBlock(const unsigned long&);
        void CopyBlock(const unsigned long& iBlock, zdouble *block);
        unsigned long GetBlockCount(void) const;
        void SetZeroes(void);
        zdouble DotProd(const SolVector&) const;
        void InitValues(const vector<zdouble>& x0);
        vector<zdouble> ResNorm(void) const;
        void AddVector(const SolVector &vec, SolVector &result, const zdouble scale = 1) const;
        void AddVector(const SolVector &vec, const zdouble scale = 1) const;
        void SubtractVector(const SolVector &vec, SolVector &result, const zdouble scale = 1) const;
        void SubtractVector(const SolVector &vec, const zdouble scale = 1);
        zdouble ComputeMagnitude(void) const;
        zdouble SquaredNorm(void) const;
        zdouble Norm(void) const;

        void WriteValues(string filename) const;

    private:
        zdouble *val = nullptr;
        unsigned short nVar;
        unsigned short nBlock;
        unsigned long len;
        void Init(void);
};
