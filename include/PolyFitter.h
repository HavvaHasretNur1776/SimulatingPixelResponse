#pragma once

#include <optional>
#include <vector>

#include "myCholeskyDecomp.h"
#include "span.h"
#include "PackedSymmetricView.h"
#include "PolyEvaluator.h"
#include "FitResult.h"

/** @brief polynomial fits
 *
 * Fits polynomials $y(dx) = \sum_{i = 1}^N c_i dx^i$ to data given data
 * (where $dx = x - x_0$ is a shift in x used to decorrelate the fit
 * parameters somewhat). The shift $x_0$ should be somewhere in center of the
 * distribution (for a linear fit, the optimal choice is
 * $x_0 = \frac{\sum_{i=1}^{N}\frac{x_i}{\sigma_{y,i}^2}}
 *             {\sum_{i=1}^{N}\frac{1}{\sigma_{y,i}^2}}$,
 * where the index $i$ runs over the measurements).
 */
template <class FLT>
class PolyFitter {
    std::size_t m_nmeas = 0;
    std::size_t m_n;
    FLT m_x0;
    std::vector<FLT> m_storage; // allocate only once for terms, rhs, mat
    
public:
    template <class FLT2>
    using Vector = cppcompat::span<FLT2>;
    template <class FLT2>
    using Matrix = PackedSymmetricView<FLT2>;

    /// constructor, support polynomials up to maxdegree, shift by x0
    PolyFitter(std::size_t maxdegree, FLT x0);

    /// clear in preparation of fitting a new data set
    void clear();
    /// return max. degree of polynomials to be used
    std::size_t maxdegree() const { return m_n - 1; }
    /// return number of measurements included this far
    std::size_t nmeas() const { return m_nmeas; }
    /// add a single data point
    void addMeas(FLT x, FLT y, FLT sigma_y);
    /// add several measurements given by iterators into corresponding ranges
    template <class IT, class JT, class KT>
    void addMeas(IT xsfirst, IT xslast, JT ysfirst, KT sigma_ysfirst)
    {
        for (; xslast != xsfirst; ++xsfirst, ++ysfirst, ++sigma_ysfirst)
            addMeas(*xsfirst, *ysfirst, *sigma_ysfirst);
    }
    /// add several measurements given by spans
    template <class FLT2>
    void addMeas(cppcompat::span<FLT2> xs, cppcompat::span<FLT2> ys,
                 cppcompat::span<FLT2> sigma_ys)
    {
        assert(xs.size() == ys.size() && xs.size() == sigma_ys.size());
        addMeas(xs.begin(), xs.end(), ys.begin(), sigma_ys.begin());
    }

    using result_type = FitResult<FLT, PolyEvaluator<FLT>>;
    std::optional<result_type> fit(std::size_t deg) const;

    // these are for debugging purposes
    Vector<const FLT> terms() const
    {
        return Vector<const FLT>{m_storage.data(), m_n};
    }
    Vector<const FLT> rhs() const
    {
        return Vector<const FLT>{m_storage.data() + m_n, m_n};
    }
    Matrix<const FLT> mat() const
    {
        return Matrix<const FLT>{m_storage.data() + 2 * m_n, m_n};
    }

private:
    Vector<FLT> _terms() { return Vector<FLT>{m_storage.data(), m_n}; }
    Vector<FLT> _rhs() { return Vector<FLT>{m_storage.data() + m_n, m_n}; }
    Matrix<FLT> _mat() { return Matrix<FLT>{m_storage.data() + 2 * m_n, m_n}; }
};

template <class FLT>
PolyFitter<FLT>::PolyFitter(std::size_t maxdegree, FLT x0)
        : m_n(1 + maxdegree), m_x0(x0),
          m_storage((m_n * (m_n + 1)) / 2 + 2 * m_n, FLT(0))
{}

template <class FLT>
void PolyFitter<FLT>::clear()
{
    m_nmeas = 0;
    std::fill(m_storage.begin() + m_n, m_storage.end(), FLT(0));
}

template <class FLT>
void PolyFitter<FLT>::addMeas(FLT x, FLT y, FLT sigma_y)
{
    auto terms_ = _terms();
    {
        auto tmp = FLT(1) / sigma_y;
        const auto dx = x - m_x0;
        for (auto& el : terms_) {
            el = tmp;
            tmp *= dx;
        }
    }
    const auto n = m_n;
    {
        y *= terms_[0];
        auto rhs_ = _rhs();
        auto mat_ = _mat();
        for (std::size_t i = 0; n != i; ++i) {
            const auto ti = terms_[i];
            for (std::size_t j = 0; j <= i; ++j) mat_[i][j] += ti * terms_[j];
            rhs_[i] += ti * y;
        }
    }
    ++m_nmeas;
}

template <class FLT>
std::optional<typename PolyFitter<FLT>::result_type>
PolyFitter<FLT>::fit(std::size_t deg) const
{
    // quick check if enough measurements for degree
    const auto n = 1 + deg, nmeas_ = nmeas();
    std::cout<<"m_n (matrixSize): "<<m_n<<" n (1+degree): "<<n<<" nmeas: "<<nmeas_<<std::endl;
    if (nmeas_ <= deg || n > m_n) return {};
    //std::cout<<__PRETTY_FUNCTION__ << " line " << __LINE__<<std::endl; 
    my::ROOT::Math::CholeskyDecompGenDim<FLT> decomp(n, mat());
    //std::cout<<__PRETTY_FUNCTION__ << " line " << __LINE__<<std::endl;
    if (!decomp) return {};
    //std::cout<<__PRETTY_FUNCTION__ << " line " << __LINE__<<std::endl;
    return result_type{n, nmeas_, std::move(decomp), rhs(),
                       PolyEvaluator<FLT>{m_x0}};
}

// vim: sw=4:tw=78:et:ft=cpp
