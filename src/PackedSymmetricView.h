#pragma once

#include <algorithm>
#include <cassert>

template <class FLT>
class PackedSymmetricView {
private:
    FLT* m_first = nullptr;
    std::size_t m_n = 0;

public:
    PackedSymmetricView(FLT* first, std::size_t n) : m_first(first), m_n(n) {}

    std::size_t n() const { return m_n; }

    FLT& operator()(std::size_t i, std::size_t j)
    {
        assert(i < m_n);
        assert(j < m_n);
        if (j > i) std::swap(i, j);
        return m_first[(i * (i + 1)) / 2 + j];
    }

    const FLT& operator()(std::size_t i, std::size_t j) const
    {
        assert(i < m_n);
        assert(j < m_n);
        if (j > i) std::swap(i, j);
        return m_first[(i * (i + 1)) / 2 + j];
    }

private:
    template <class PARENT> friend class rowproxy;

    template <class PARENT>
    class rowproxy {
    private:
        PARENT& m_parent;
        std::size_t m_i;

    public:
        rowproxy(PARENT& parent, std::size_t i) : m_parent(parent), m_i(i) {}

        auto& operator[](std::size_t j)
        {
            assert(j < m_parent.m_n);
            std::size_t i = m_i;
            if (j > i) std::swap(i, j);
            return m_parent.m_first[(i * (i + 1)) / 2 + j];
        }
        const auto& operator[](std::size_t j) const
        {
            assert(j < m_parent.m_n);
            std::size_t i = m_i;
            if (j > i) std::swap(i, j);
            return m_parent.m_first[(i * (i + 1)) / 2 + j];
        }
    };

public:
    rowproxy<PackedSymmetricView> operator[](std::size_t i)
    {
        assert(i < m_n);
        return {*this, i};
    }
    rowproxy<const PackedSymmetricView> operator[](std::size_t i) const
    {
        assert(i < m_n);
        return {*this, i};
    }
};

// vim: sw=4:tw=78:et:ft=cpp
