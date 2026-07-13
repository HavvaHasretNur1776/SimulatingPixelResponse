#pragma once

#include <numeric>

/// evaluate a polynomial given its coefficients
template <class FLT>
class PolyEvaluator {
    FLT m_x0;

public:
    /// constructor
    constexpr explicit PolyEvaluator(FLT x0 = 0.) noexcept : m_x0(x0) {}
    /// evaluate polynomial given x and coefficients
    template <class COEFFS>
    constexpr FLT operator()(FLT x, const COEFFS& coeffs) noexcept
    {
        const auto dx = x - m_x0;
        return std::accumulate(
                coeffs.rbegin(), coeffs.rend(), FLT(0),
                [dx](auto sum, auto coeff) { return coeff + dx * sum; });
    }
};

// vim: sw=4:tw=78:et:ft=cpp
