/**
 * @author Tim Luchterhand
 * @date 21.03.23.
 */

#ifndef TEMPO_MATRIX_HPP
#define TEMPO_MATRIX_HPP

#include <cmath>
#include <functional>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <Iterators.hpp>
#include <nlohmann/json.hpp>

#include "util/traits.hpp"

namespace tempo {

    /**
     * Matrix storage layout types
     */
    enum class Layout {
        RowMajor, ColMajor
    };

    namespace detail {
        template<typename Iterator>
        class StrideIterator : public iterators::impl::SynthesizedOperators<StrideIterator<Iterator>> {
            Iterator it;
            Iterator::difference_type stride;
        public:
            using value_type = typename Iterator::value_type;
            using reference = typename Iterator::reference;
            using pointer = typename Iterator::pointer;
            using difference_type = typename Iterator::difference_type;
            using iterator_category = typename Iterator::iterator_category;
            using iterator_concept = typename Iterator::iterator_concept;
            using iterator_type = typename Iterator::iterator_type;
            using iterators::impl::SynthesizedOperators<StrideIterator>::operator++;
            using iterators::impl::SynthesizedOperators<StrideIterator>::operator--;

            constexpr StrideIterator() noexcept = default;

            constexpr StrideIterator(Iterator it, difference_type stride) noexcept: it(it), stride(stride) {}

            template<typename I>
            constexpr StrideIterator(const StrideIterator<I> &other) noexcept :
                    it(other.getIterator()), stride(other.getStride()) {}

            Iterator getIterator() const noexcept {
                return it;
            }

            auto getStride() const noexcept {
                return stride;
            }

            constexpr StrideIterator &operator++() noexcept {
                it += stride;
                return *this;
            }

            constexpr StrideIterator &operator--() noexcept {
                it -= stride;
                return *this;
            }

            constexpr StrideIterator &operator+=(difference_type n) noexcept {
                it += stride * n;
                return *this;
            }

            constexpr StrideIterator &operator-=(difference_type n) noexcept {
                it -= stride * n;
                return *this;
            }

            constexpr difference_type operator-(const StrideIterator &other) const noexcept {
                return static_cast<difference_type>((it - other.it) / other.stride);
            }

            constexpr bool operator==(const StrideIterator &other) const noexcept {
                return it == other.it;
            }

            constexpr reference operator*() const noexcept {
                return *it;
            }

            constexpr reference operator*() noexcept {
                return *it;
            }

            constexpr auto operator->() const noexcept {
                return it.operator->();
            }

            constexpr auto operator->() noexcept {
                return it.operator->();
            }

            constexpr bool operator<(const StrideIterator &other) const noexcept {
                return it < other.it;
            }

            constexpr bool operator>(const StrideIterator &other) const noexcept {
                return it > other.it;
            }
        };
    }


    /**
     * @brief Very simple implementation of a matrix using a single vector
     * @tparam T element type
     */
    template<typename T>
    class Matrix {

        [[nodiscard]] constexpr std::size_t index(std::size_t i, std::size_t j) const noexcept {
            return i * (layout == Layout::RowMajor ? nCols : 1) + j * (layout == Layout::RowMajor ? 1 : nRows);
        }

        [[nodiscard]] unsigned printLength() const noexcept {
            unsigned max = 0;
            for_each([&max](const auto &val) {
                if constexpr (std::integral<T>) {
                    auto decimals = static_cast<unsigned>(std::log(static_cast<float>(std::abs(val)) + 1)
                            / std::log(10)) + (val < 0);
                    max = std::max(max, decimals);
                } else {
                    std::stringstream ss;
                    ss << val;
                    max = std::max(max, static_cast<unsigned>(ss.str().length()));
                }
            });

            return max + 2;
        }

        constexpr auto rowIterators(std::size_t r) noexcept {
            std::ptrdiff_t stride = layout == Layout::RowMajor ? 1 : nRows;
            detail::StrideIterator begin(data.begin() + index(r, 0), stride);
            detail::StrideIterator end(data.begin() + index(r, nCols), stride);
            return std::make_pair(begin, end);
        }

        constexpr auto colIterators(std::size_t c) noexcept {
            std::ptrdiff_t stride = layout == Layout::RowMajor ? nCols : 1;
            detail::StrideIterator begin(data.begin() + index(0, c), stride);
            detail::StrideIterator end(data.begin() + index(nRows, c), stride);
            return std::make_pair(begin, end);
        }

        std::vector<T> data;
        std::size_t nRows{}, nCols{};
        Layout layout{};

    public:
        using reference = typename std::vector<T>::reference;
        using const_reference = typename std::vector<T>::const_reference;
        using value_type = typename std::vector<T>::value_type;

        /**
         * Default CTor. Constructs an empty Matrix
         */
        constexpr Matrix() = default;

        constexpr void swap(Matrix &other) noexcept {
            using std::swap;
            swap(data, other.data);
            swap(nRows, other.nRows);
            swap(nCols, other.nCols);
            swap(layout, other.layout);
        }

        constexpr Matrix(const Matrix &other) = default;

        constexpr Matrix(Matrix &&other) noexcept: Matrix() {
            swap(other);
        }

        constexpr Matrix &operator=(Matrix other) noexcept {
            swap(other);
            return *this;
        }

        constexpr ~Matrix() = default;


        /**
         * CTor. Initializes a matrix with a given value
         * @param nRows number of rows
         * @param nCols  number of columns
         * @param initValue initial value (optional, if non is specified, default constructs elements)
         * @param layout storage layout (default is row major)
         */
        constexpr Matrix(std::size_t nRows, std::size_t nCols, const T &initValue = T(),
                         Layout layout = Layout::RowMajor) : data(nCols * nRows, initValue), nRows(nRows), nCols(nCols),
                                                             layout(layout) {}

        /**
         * Ctor. Initializes a matrix with the values from the given range
         * @tparam R type of range
         * @param nRows number of rows
         * @param nCols  number of columns
         * @param values range of values to
         * @param layout storage layout (default is row major)
         */
        template<std::ranges::range R>
        constexpr Matrix(std::size_t nRows, std::size_t nCols, const R &values, Layout layout = Layout::RowMajor):
                data(std::ranges::begin(values), std::ranges::end(values)), nRows(nRows), nCols(nCols),
                layout(layout) {

            if (data.size() != nRows * nCols) {
                throw std::runtime_error("wrong number of values for matrix initialization.");
            }
        }

        /**
         * Ctor. Initializes a matrix with the values from the given initializer list
         * @param nRows number of rows
         * @param nCols  number of columns
         * @param values range of values to
         * @param layout storage layout (default is row major)
         */
        constexpr Matrix(std::size_t nRows, std::size_t nCols, std::initializer_list<T> values,
                         Layout layout = Layout::RowMajor) : data(values.begin(), values.end()), nRows(nRows),
                                                             nCols(nCols), layout(layout) {
            if (data.size() != nRows * nCols) {
                throw std::runtime_error("wrong number of values for matrix initialization.");
            }
        }

        /**
         * @tparam F matrix like type that supports function call operator
         * @param nRows number of rows
         * @param nCols  number of columns
         * @param matrixLike functor object with init values
         * @param initValue initial value (optional, if non is specified, default constructs elements)
         * @param layout storage layout (default is row major)
         */
        template<concepts::callable_r<T, std::size_t, std::size_t> F>
        constexpr Matrix(std::size_t nRows, std::size_t nCols, F &&matrixLike, Layout layout = Layout::RowMajor):
                data(nRows * nCols), nRows(nRows), nCols(nCols), layout(layout) {
            if (layout == Layout::RowMajor) {
                for (std::size_t i = 0; i < numRows(); ++i) {
                    for (std::size_t j = 0; j < numColumns(); ++j) {
                        operator()(i, j) = std::forward<F>(matrixLike)(i, j);
                    }
                }
            } else {
                for (std::size_t j = 0; j < numColumns(); ++j) {
                    for (std::size_t i = 0; i < numRows(); ++i) {
                        operator()(i, j) = std::forward<F>(matrixLike)(i, j);
                    }
                }
            }
        }

        /**
         * Resize the matrix. Does not change the stored items. Depending on the new size, old items might be lost
         * @param numRows new number of rows
         * @param numCols new number of columns
         * @param value value to insert in case of insertion
         */
        constexpr void resize(std::size_t numRows, std::size_t numCols, const T& value = T{}) {
            data.resize(numRows * numCols, value);
            nRows = numRows;
            nCols = numCols;
        }

        /**
         * Fill the matrix with the given value
         * @param value value to insert at every position
         */
        constexpr void fill(const T &value) {
            for_each([&value](auto &val) { val = value; });
        }

        /**
         * Set the size of the matrix and insert the given value at every position
         * @param numRows number of rows
         * @param numCols number of columns
         * @param value value to insert at every position
         * @note in contrast to resize this method overwrites existing elements
         */
        constexpr void fill(std::size_t numRows, std::size_t numCols, const T &value) {
            std::size_t newSize = numRows * numCols;
            nRows = numRows;
            nCols = numCols;
            if (newSize > data.size()) {
                data = std::vector<T>(newSize, value);
            } else {
                data.resize(newSize);
                fill(value);
            }
        }

        /**
         * Change the storage layout
         * @param newLayout new layout
         */
        constexpr void changeLayout(Layout newLayout) noexcept {
            layout = newLayout;
        }

        /**
         * Element access without bounds checking
         * @param i row index
         * @param j column index
         * @return reference to value
         */
        constexpr reference operator()(std::size_t i, std::size_t j) noexcept {
            return data[index(i, j)];
        }

        /**
         * @copydoc operator()()
         */
        constexpr const_reference operator()(std::size_t i, std::size_t j) const noexcept {
            return data[index(i, j)];
        }

        /**
         * Element access with bounds checking
         * @param i row index
         * @param j column index
         * @return reference to value
         * @throws std::out_of_range if index out of range
         */
        constexpr reference at(std::size_t i, std::size_t j) {
            if (i >= nRows || j >= nCols) {
                throw std::out_of_range(std::string(
                        "access element at (" + std::to_string(i) + ", " + std::to_string(j) +
                        ") but matrix has dimensions (" + std::to_string(nRows) + ", " + std::to_string(nCols) + ")"));
            }

            return (*this)(i, j);
        }

        /**
         * \copydoc at
         */
        constexpr const_reference at(std::size_t i, std::size_t j) const {
            return const_cast<Matrix *>(this)->at(i, j);
        }

        /**
         * Get the number of rows
         * @return
         */
        [[nodiscard]] constexpr std::size_t numRows() const noexcept {
            return nRows;
        }

        /**
         * Get the number of columns
         * @return
         */
        [[nodiscard]] constexpr std::size_t numColumns() const noexcept {
            return nCols;
        }

        /**
         * get the current storage layout
         * @return
         */
        [[nodiscard]] constexpr Layout storageLayout() const noexcept {
            return layout;
        }

        /**
         * Applies a given functor to all elements of the matrix
         * @param functor callable that takes a mutable reference to a matrix element
         */
        constexpr void for_each(std::invocable<T &> auto &&functor) {
            for(auto &&val : data) {
                std::forward<decltype(functor)>(functor)(std::forward<decltype(val)>(val));
            }
        }

        /**
         * @copydoc for_each
         */
        constexpr void for_each(std::invocable<T &> auto &&functor) const {
            for(auto &&val : data) {
                std::forward<decltype(functor)>(functor)(std::forward<decltype(val)>(val));
            }
        }

        /**
         * reference to underlying data structure
         * @return
         */
        constexpr const auto &rawData() const noexcept {
            return data;
        }

        /**
         * @copydoc rawData
         */
        constexpr auto &rawData() noexcept {
            return data;
        }

        template<typename U = T> requires(concepts::printable<U>)
        friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
            auto printLen = matrix.printLength();
            os << std::setprecision(2);
            for (std::size_t i = 0; i < matrix.numRows(); ++i) {
                for (std::size_t j = 0; j < matrix.numColumns(); ++j) {
                    os << std::setw(printLen) << matrix(i, j);
                }

                os << std::setw(0) << std::endl;
            }

            return os;
        }

        /**
         * equality comparison
         * @param other
         * @return returns true if the matrices have the same dimensions and the same elements (independent of storage
         * layout), else false
         */
        constexpr bool operator==(const Matrix &other) const noexcept {
            if (numRows() != other.numRows() or numColumns() != other.numColumns()) {
                return false;
            }

            for (auto r = 0ul; r < numRows(); ++r) {
                for (auto c = 0ul; c < numColumns(); ++c) {
                    if ((*this)(r, c) != other(r, c)) {
                        return false;
                    }
                }
            }

            return true;
        }

        /**
         * Row view of the matrix without bounds checking
         * @param r row index
         * @return row view of the specified row
         */
        constexpr auto row_unsafe(std::size_t r) noexcept {
            auto [begin, end] = rowIterators(r);
            return std::ranges::subrange(begin, end);
        }

        /**
         * @copydoc row_unsafe()
         */
        constexpr auto row_unsafe(std::size_t r) const noexcept {
            using It = detail::StrideIterator<typename std::vector<T>::const_iterator>;
            auto [begin, end] = const_cast<Matrix*>(this)->rowIterators(r);
            return std::ranges::subrange(It(begin), It(end));
        }

        /**
         * Column view of the matrix without bounds checking
         * @param c column index
         * @return column view of the specified column
         */
        constexpr auto column_unsafe(std::size_t c) noexcept {
            auto [begin, end] = colIterators(c);
            return std::ranges::subrange(begin, end);
        }

        /**
         * @copydoc column_unsafe()
         */
        constexpr auto column_unsafe(std::size_t c) const noexcept {
            using It = detail::StrideIterator<typename std::vector<T>::const_iterator>;
            auto [begin, end] = const_cast<Matrix*>(this)->colIterators(c);
            return std::ranges::subrange(It(begin), It(end));
        }

        /**
         * Row view of the matrix
         * @param r row index
         * @return row view of the specified row
         * @throws std::out_of_range if row index is out of bounds
         */
        constexpr auto row(std::size_t r) {
            if (r >= nRows) {
                throw std::out_of_range(std::string(
                        "access row " + std::to_string(r) +
                        " but matrix has dimensions (" + std::to_string(nRows) + ", " + std::to_string(nCols) + ")"));
            }

            return row_unsafe(r);
        }

        /**
         * @copydoc row()
         */
        constexpr auto row(std::size_t r) const {
            if (r >= nRows) {
                throw std::out_of_range(std::string(
                        "access row " + std::to_string(r) +
                        " but matrix has dimensions (" + std::to_string(nRows) + ", " + std::to_string(nCols) + ")"));
            }

            return row_unsafe(r);
        }

        /**
         * Column view of the matrix
         * @param c column index
         * @return column view of the specified column
         * @throws std::out_of_range if column index is out of bounds
         */
        constexpr auto column(std::size_t c) {
            if (c >= nCols) {
                throw std::out_of_range(std::string(
                        "access col " + std::to_string(c) +
                        " but matrix has dimensions (" + std::to_string(nRows) + ", " + std::to_string(nCols) + ")"));
            }

            return column_unsafe(c);
        }

        /**
         * @copydoc column()
         */
        constexpr auto column(std::size_t c) const {
            if (c >= nCols) {
                throw std::out_of_range(std::string(
                        "access col " + std::to_string(c) +
                        " but matrix has dimensions (" + std::to_string(nRows) + ", " + std::to_string(nCols) + ")"));
            }

            return column_unsafe(c);
        }
    };
}

namespace nlohmann {
    template<typename T>
    struct adl_serializer<tempo::Matrix<T>> {
        static void to_json(json &j, const tempo::Matrix<T> &matrix) {
            using Layout = tempo::Layout;
            j["data"] = matrix.rawData();
            j["nRows"] = matrix.numRows();
            j["nCols"] = matrix.numColumns();
            j["layout"] = matrix.storageLayout() == Layout::RowMajor ? "RowMajor" : "ColMajor";
        }

        static void from_json(const json &j, tempo::Matrix<T> &matrix) {
            matrix = tempo::Matrix<T>(j.at("nRows"), j.at("nCols"), j.at("data"),
                                      j.at("layout").get<std::string>() == "RowMajor" ? tempo::Layout::RowMajor
                                                                                      : tempo::Layout::ColMajor);
        }
    };
}
#endif //TEMPO_MATRIX_HPP
