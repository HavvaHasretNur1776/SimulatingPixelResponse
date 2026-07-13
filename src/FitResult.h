#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "myCholeskyDecomp.h" // "Math/CholeskyDecomp.h"

#include "span.h"
#include "PackedSymmetricView.h"

template <class FLT, class EVALUATOR>
class FitResult {
public:
    using Decomp = my::ROOT::Math::CholeskyDecompGenDim<FLT>;

private:
    std::size_t m_n;
    std::size_t m_nmeas;
    Decomp m_decomp;
    std::vector<FLT> m_storage;
    bool m_solved = false;
    EVALUATOR m_eval;

public:
    template <class FLT2>
    using Vector = cppcompat::span<FLT2>;
    template <class FLT2>
    using Matrix = PackedSymmetricView<FLT2>;

    template <class RHS>
    FitResult(std::size_t n, std::size_t nmeas, Decomp&& decomp,
              const RHS& rhs, EVALUATOR&& eval)
            : m_n(n), m_nmeas(nmeas), m_decomp(std::move(decomp)),
              m_eval(std::move(eval))
    {
        m_storage.reserve(n + (n * (n + 1)) / 2);
        m_storage.insert(m_storage.end(), rhs.begin(), rhs.begin() + n);
        // ... and pad with NaNs to the size of the covariance matrix (which
        // allows us to recognize when we have not inverted the matrix yet)
        m_storage.resize(n + (n * (n + 1)) / 2, std::numeric_limits<FLT>::quiet_NaN());
    }

    std::size_t nparams() const noexcept { return m_n; }
    std::size_t nmeas() const noexcept { return m_nmeas; }
    std::size_t ndf() const noexcept { return nmeas() - nparams(); }

    operator bool() const noexcept { return !m_storage.empty(); }
    bool operator!() const noexcept { return m_storage.empty(); }

    Vector<const FLT> params() noexcept
    {
        if (!m_solved) {
            // not solved yet, do so now
            Vector<FLT> rhs{m_storage.data(), m_n};
            m_decomp.Solve(rhs);
            m_solved = true;
        }
        return Vector<const FLT>{m_storage.data(), m_n};
    }

    Matrix<const FLT> cov() noexcept
    {
        if (std::isnan(m_storage[m_n])) {
            // not inverted yet, do so now
            Matrix<FLT> mat{m_storage.data() + m_n, m_n};
            m_decomp.Invert(mat);
        }
        return Matrix<const FLT>{m_storage.data() + m_n, m_n};
    }
    ///Get Chi2 Matrix
    Matrix<const FLT> matrix() noexcept
    {
       //Returning elements before inverting
            //if solved Matrix would be cov() 
        return Matrix<const FLT>{m_storage.data() + m_n, m_n};
    }
    ////Get RHS vector
    Vector<const FLT> rhs() noexcept
    {
        //Returning elements before solving
            //if solved rhs would be params() 
        return Vector<const FLT>{m_storage.data(), m_n};
    }

    FLT chi2(FLT x, FLT y, FLT sigma_y)
    {
        //const auto dy = (*y - m_eval(*x, params())) / *sigma_y;
        const auto dy = (y - m_eval(x, params())) / sigma_y;

        return dy * dy;
    }

    template <class IT, class JT, class KT>
    FLT chi2(IT xsfirst, IT xslast, JT ysfirst, KT sigma_ysfirst)
    {
        const auto params_ = params();
        FLT retVal = 0;
        for (; xslast != xsfirst; ++xsfirst, ++ysfirst, ++sigma_ysfirst) {
            const auto dy =
                    (*ysfirst - m_eval(*xsfirst, params_)) / *sigma_ysfirst;
            retVal += dy * dy;
        }
        return retVal;
    }

    template <class FLT2>
    FLT chi2(cppcompat::span<FLT2> xs, cppcompat::span<FLT2> ys,
             cppcompat::span<FLT2> sigma_ys)
    {
        assert(xs.size() == ys.size() && xs.size() == sigma_ys.size());
        return chi2(xs.begin(), xs.end(), ys.begin(), sigma_ys.begin());
    }

////Get fitted vector
    FLT fittedCoordinate(FLT x)
    //Vector<const FLT> fittedCoordinate()
    {
        const auto params_ = params();
        return m_eval(x, params_);
    } 


};

// vim: sw=4:tw=78:et:ft=cpp
