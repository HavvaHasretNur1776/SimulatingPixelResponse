#pragma once

#if __cplusplus >= 202002L
#include <span>
namespace cppcompat {
    // we have C++20 or later, so use the real deal
    using std::dynamic_extent;
    using std::span;
} // namespace cppcompat

#else // __cplusplus >= 202002L

#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <deque>
#include <type_traits>

namespace cppcompat {
    /// used to indicate that length of span may only be known at run time
    constexpr std::size_t dynamic_extent =
            std::numeric_limits<std::size_t>::max();

    // forward declaration
    template <class T, std::size_t EXTENT>
    class span;

    /// namespace for implementation details
    namespace impl {
        // get size of containers in a uniform way
        template <class CONTAINER>
        constexpr std::size_t size(const CONTAINER& c) noexcept
        {
            return c.size();
        }
        template <class T, std::size_t SZ>
        constexpr std::size_t size(const T (&arr)[SZ]) noexcept
        {
            static_cast<void>(arr);
            return SZ;
        }
        template <class T, std::size_t SZ>
        constexpr std::size_t size(const volatile T (&arr)[SZ]) noexcept
        {
            static_cast<void>(arr);
            return SZ;
        }

        // get contiguous data from container in uniform way
        template <class CONTAINER>
        constexpr auto data(CONTAINER& c) noexcept -> decltype(c.data())
        {
            return c.data();
        }
        template <class T, std::size_t SZ>
        constexpr T* data(T (&arr)[SZ]) noexcept
        {
            return arr;
        }

        // forward declaration
        template <class ELEMENT_TYPE, std::size_t EXTENT, class CONTAINER>
        struct is_compat_container;

        // helper for span's size method
        template <std::size_t SZ>
        struct size_helper {
            constexpr static std::size_t
            size(const std::size_t& /* unused */) noexcept
            {
                return SZ;
            }
        };
        template <>
        struct size_helper<dynamic_extent> {
            constexpr static const std::size_t&
            size(const std::size_t& sz) noexcept
            {
                return sz;
            }
        };
    } // namespace impl

    /** @brief Represent span of Ts allocated contiguously in memory.
     *
     * A span represents a span of Ts which is allocated contiguously in
     * memory. It can either be a span whose size is known at compile-time
     * (EXTENT != dynamic_extent) or a span whose size is only known at
     * run-time (EXTENT == dynamic_extent).
     *
     * It is a near drop-in replacement for C++20's span (up to things C++14
     * cannot do, and bugs and limitations due to the author's [in]abilities).
     * Since C++14 lacks class template argument deduction guides, this
     * implementation offers a make_span function that helps somewhat to work
     * around that lack.
     *
     * A span can be used to erase the type of a container/array of elements
     * for containers that are laid out contiguously in memory. Here's a brief
     * example:
     *
     * @code
     * template <class T>
     * void square(span<T> dst, const span<const T> src)
     * {
     *     assert(dst.size() == src.size());
     *     auto it = std::begin(dst);
     *     for (const auto& el: src) *it++ = el * el;
     * }
     *
     * void useSquare()
     * {
     *     const std::array<int, 4> src{{1, 2, 3, 4}};
     *     std::vector<int> dst1(4);
     *     std::array<int, 4> dst2;
     *     int dst3[4];
     *     square(make_span(dst1), make_span(src));
     *     square(make_span(dst2), make_span(src));
     *     square(make_span(dst3), make_span(src));
     * }
     * @endcode
     */
    template <class T, std::size_t EXTENT = dynamic_extent>
    class span {
    public:
        using element_type = T;
        using value_type = typename std::remove_cv<element_type>::type;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using const_pointer = const T*;
        using reference = T&;
        using const_reference = const T&;
        using iterator = pointer;
        using reverse_iterator = std::reverse_iterator<iterator>;
        static constexpr size_type extent = EXTENT;

    private:
        size_type m_sz = 0;
        pointer m_ptr = nullptr;

    public:
        constexpr span() noexcept = default;
        constexpr span(const span&) noexcept = default;
        constexpr span& operator=(const span&) noexcept = default;
        ~span() noexcept = default;

        constexpr explicit span(pointer first, size_type size) noexcept
                : m_sz((assert(EXTENT == dynamic_extent || EXTENT == size),
                        size)),
                  m_ptr(first)
        {}
        constexpr explicit span(pointer first, pointer last) noexcept
                : span(first, std::distance(first, last))
        {}

        // this constructor takes care of allowing conversions of spans if you
        // create a span<const T> from a span<T>, and if you change from/to
        // dynamic_extent
        template <class TT, std::size_t SZ,
                  class = typename std::enable_if<
                          std::is_convertible<TT*, element_type*>::value &&
                          (SZ == EXTENT || SZ == dynamic_extent ||
                           EXTENT == dynamic_extent)>::type>
        constexpr span(const span<TT, SZ>& sp) noexcept
                : span(sp.begin(), sp.size())
        {}

        // takes care of construction from a "spannable" container like array
        // or vector, or from a compatible span...
        template <class CONTAINER,
                  class = typename std::enable_if<impl::is_compat_container<
                          element_type, EXTENT, CONTAINER>::value>::type>
        constexpr explicit span(CONTAINER& c) noexcept
                : span(cppcompat::impl::data(c), cppcompat::impl::size(c))
        {}

        constexpr size_type size() const noexcept
        {
            return impl::size_helper<EXTENT>::size(m_sz);
        }
        constexpr size_type size_bytes() const noexcept
        {
            return size() * sizeof(value_type);
        }
        constexpr bool empty() const noexcept { return 0 == size(); }

        constexpr pointer data() const noexcept { return m_ptr; }

        constexpr iterator begin() const noexcept { return m_ptr; }
        constexpr iterator end() const noexcept { return m_ptr + size(); }

        constexpr reverse_iterator rbegin() const noexcept
        {
            return reverse_iterator(end());
        }
        constexpr reverse_iterator rend() const noexcept
        {
            return reverse_iterator(begin());
        }

        constexpr reference front() const noexcept
        {
            assert(!empty());
            return *begin();
        }
        constexpr reference back() const noexcept
        {
            assert(!empty());
            return *(begin() + (size() - 1));
        }

        constexpr reference operator[](size_type idx) const noexcept
        {
            assert(idx < size());
            return *(begin() + idx);
        }

        template <size_type SZ>
        constexpr span<element_type, SZ> first() const noexcept
        {
            assert(SZ <= size());
            return span<element_type, SZ>{begin(), SZ};
        }
        template <size_type SZ>
        constexpr span<element_type, SZ> last() const noexcept
        {
            assert(SZ <= size());
            return span<element_type, SZ>{begin() + (size() - SZ), SZ};
        }
        template <size_type OFFSET, size_type SZ>
        constexpr span<element_type, SZ> subspan() const noexcept
        {
            assert(SZ <= OFFSET + SZ && OFFSET <= OFFSET + SZ);
            assert(OFFSET <= size());
            assert(OFFSET + SZ <= size());
            return span<element_type, SZ>{begin() + OFFSET, SZ};
        }

        constexpr span<element_type> first(const size_type sz) const noexcept
        {
            assert(sz <= size());
            return span<element_type>{begin(), sz};
        }
        constexpr span<element_type> last(const size_type sz) const noexcept
        {
            assert(sz <= size());
            return span<element_type>{begin() + (size() - sz), sz};
        }
        constexpr span<element_type>
        subspan(const size_type offset,
                size_type sz = dynamic_extent) const noexcept
        {
            assert(offset <= size());
            if (dynamic_extent == sz) sz = size() - offset;
            assert(sz <= offset + sz && offset <= offset + sz);
            assert(offset + sz <= size());
            return span<element_type>{begin() + offset, sz};
        }
    };

    static_assert(std::is_trivially_copyable<span<int>>::value,
                  "span implementation error");

    namespace impl {
        // helpers
        template <class ELEMENT_TYPE, class CONTAINER,
                  class PTR_TYPE = decltype(cppcompat::impl::data(
                          std::declval<CONTAINER&>()))>
        constexpr std::is_convertible<PTR_TYPE, ELEMENT_TYPE*>
        _has_compat_data(int /* unused */) noexcept;
        template <class ELEMENT_TYPE, class CONTAINER>
        constexpr std::false_type
        _has_compat_data(long /* unused */) noexcept;
        /// check if container has data that is compatible with ELEMENT_TYPE
        template <class ELEMENT_TYPE, class CONTAINER>
        constexpr bool has_compat_data() noexcept
        {
            return decltype(_has_compat_data<ELEMENT_TYPE, CONTAINER>(
                    0))::value;
        }

        static_assert(!has_compat_data<int, int>(), "has_compat_data broken");
        static_assert(has_compat_data<int, int (&)[3]>(),
                      "has_compat_data broken");
        static_assert(has_compat_data<int, std::array<int, 1>&>(),
                      "has_compat_data broken");
        static_assert(has_compat_data<const int, std::array<int, 1>&>(),
                      "has_compat_data broken");
        static_assert(!has_compat_data<int, const std::array<int, 1>&>(),
                      "has_compat_data broken");
        static_assert(has_compat_data<int, span<int, 1>>(),
                      "has_compat_data broken");
        static_assert(has_compat_data<const int, span<int, 1>>(),
                      "has_compat_data broken");

        /// compile-time size of container or dynamic_extent if unavailable
        template <class CONTAINER>
        struct _size : std::integral_constant<std::size_t, dynamic_extent> {};
        template <class T, std::size_t EXTENT>
        struct _size<span<T, EXTENT>>
                : std::integral_constant<std::size_t, EXTENT> {};
        template <class T, std::size_t EXTENT>
        struct _size<std::array<T, EXTENT>>
                : std::integral_constant<std::size_t, EXTENT> {};
        template <class T, std::size_t EXTENT>
        struct _size<T (&)[EXTENT]>
                : std::integral_constant<std::size_t, EXTENT> {};

        /// is CONTAINER compatible with span<ELEMENT_TYPE, EXTENT>
        template <class ELEMENT_TYPE, std::size_t EXTENT, class CONTAINER>
        struct is_compat_container
                : std::integral_constant<
                          bool, has_compat_data<ELEMENT_TYPE, CONTAINER>() &&
                                        (dynamic_extent == EXTENT ||
                                         EXTENT <= _size<CONTAINER>::value)> {
        };

        static_assert(is_compat_container<int, 3, std::array<int, 3>>::value,
                      "is_compat_container is broken");
        static_assert(is_compat_container<int, 2, std::array<int, 3>>::value,
                      "is_compat_container is broken");
        static_assert(!is_compat_container<int, 4, std::array<int, 3>>::value,
                      "is_compat_container is broken");
        static_assert(
                !is_compat_container<long, 3, std::array<int, 3>>::value,
                "is_compat_container is broken");
        static_assert(is_compat_container<int, dynamic_extent,
                                          std::array<int, 3>>::value,
                      "is_compat_container is broken");
        static_assert(is_compat_container<int, 2, span<int, 3>>::value,
                      "is_compat_container is broken");
        static_assert(
                is_compat_container<int, 2, span<int, dynamic_extent>>::value,
                "is_compat_container is broken");
    } // namespace impl
} // namespace cppcompat

#endif // __cplusplus >= 202002L

namespace cppcompat {
    // make_span is not part of C++20, but we need it to work around C++14's
    // missing class template argument deduction guides
    template <class T, std::size_t N>
    span<T, N> make_span(T (&arr)[N]) noexcept
    {
        return span<T, N>{arr};
    }

    template <class T, std::size_t N>
    span<T, N> make_span(std::array<T, N>& arr) noexcept
    {
        return span<T, N>{arr};
    }
    template <class T, std::size_t N>
    span<const T, N> make_span(const std::array<T, N>& arr) noexcept
    {
        return span<const T, N>{arr};
    }

    template <class CONTAINER>
    span<typename CONTAINER::value_type, dynamic_extent>
    make_span(CONTAINER& c) noexcept
    {
        return span<typename CONTAINER::value_type, dynamic_extent>{c};
    }
    template <class CONTAINER>
    span<const typename CONTAINER::value_type, dynamic_extent>
    make_span(const CONTAINER& c) noexcept
    {
        return span<const typename CONTAINER::value_type, dynamic_extent>{c};
    }

    namespace impl {
        // can't do this properly in C++17, so be a tiny bit too generous and
        // check if IT is a random access iterator, even if it's not quite
        // correct...
        template <class IT, class T= typename std::iterator_traits<IT>::value_type>
        struct is_contiguous_iterator
                : public std::is_base_of<std::random_access_iterator_tag,
                                         typename std::iterator_traits<
                                                 IT>::iterator_category> {};
        // STL's deque has random access iterators, but is not contiguous, so
        // blacklist its iterators (at least the ones of deque's using
        // std::allocator<T>)
        template <class T>
        struct is_contiguous_iterator<typename std::deque<T>::iterator, T>
                : public std::false_type {};
        template <class T>
        struct is_contiguous_iterator<typename std::deque<T>::const_iterator, T>
                : public std::false_type {};
    } // namespace impl

    template <class IT,
              class = std::enable_if_t<
                      cppcompat::impl::is_contiguous_iterator<IT>::value>>
    span<std::remove_reference_t<
                 typename std::iterator_traits<IT>::reference>,
         dynamic_extent>
    make_span(IT first, IT last) noexcept
    {
        return span<std::remove_reference_t<
                            typename std::iterator_traits<IT>::reference>,
                    dynamic_extent>{&*first, &*last};
    }
} // namespace cppcompat

// vim: sw=4:tw=78:ft=cpp:et
