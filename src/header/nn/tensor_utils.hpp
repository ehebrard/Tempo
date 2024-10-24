/**
 * @author Tim Luchterhand
 * @date 20.07.23.
 */

#ifndef SCHEDCL_TENSOR_UTILS_HPP
#define SCHEDCL_TENSOR_UTILS_HPP
#include <vector>
#include <ranges>
#include <span>
#include <Iterators.hpp>
#include <filesystem>
#include "util/traits.hpp"
#include "torch_types.hpp"

/**
 * @brief namespace containing utility functions to handle tensors
 */
namespace tempo::nn::util {

    namespace impl {
        struct SliceSurrogate {
            constexpr SliceSurrogate(IndexType val) noexcept: val(val) {}
            [[nodiscard]] constexpr bool isNone() const noexcept {
                return val < 0;
            }

            constexpr operator IndexType() const noexcept {
                return val < 0 ? 0 : val;
            }

            IndexType val;
        };
    }

    inline constexpr impl::SliceSurrogate SliceHere = -1;

    /**
     * Function used for constant indexing of scalars in tensors.
     * @tparam T underlying data type of tensor
     * @param tensor tensor to index
     * @param indices index for each dimension
     * @return const reference to tensor element
     * @note this kind of indexing operates directly on the underlying data bypassing all torch mechanisms for recording
     * gradients for example. It is thus significantly faster but also less flexible and more error prone. Make sure
     * you checked your indices! Also tensor memory needs to be contiguous.
     */
    template<concepts::scalar T>
    const T&c_index(const torch::Tensor &tensor, std::initializer_list<long> indices) {
        if (not tensor.is_contiguous()) {
            throw std::runtime_error("non-contiguous tensor cannot be indexed by this function");
        }

        if (indices.size() > tensor.sizes().size()) {
            throw std::out_of_range("to many indices");
        }

        std::ptrdiff_t offset = 0;
        for (auto [stride, index] : iterators::const_zip(tensor.strides(), indices)) {
            offset += stride * index;
        }

        if (offset < 0 or offset > tensor.numel()) {
            throw std::out_of_range("out of range access to tensor storage");
        }

        return tensor.data_ptr<T>()[offset];
    }

    /**
     * index version that returns mutable reference to tensor element (see c_index)
     */
    template<concepts::scalar T>
    T &index(torch::Tensor &tensor, std::initializer_list<long> indices) {
        return const_cast<T&>(c_index<T>(tensor, indices));
    }

    /**
     * Assigns values to a 1D slice of a tensor
     * @tparam T underlying data type of tensor
     * @param tensor tensor to modify
     * @param indices index for each dimension. A negative value indicates the dimension along which to slice. This
     * supports slicing only along ONE SINGLE DIMENSION!
     * @param values values to insert into the slice
     */
    template<concepts::scalar T>
    auto sliceAssign(torch::Tensor &tensor, std::initializer_list<impl::SliceSurrogate> indices,
                     std::initializer_list<T> values) -> torch::Tensor & {
        if (not tensor.is_contiguous()) {
            throw std::runtime_error("non-contiguous tensor cannot be indexed by this function");
        }

        if (indices.size() > tensor.sizes().size()) {
            throw std::out_of_range("to many indices");
        }

        std::ptrdiff_t baseOffset = 0;
        std::ptrdiff_t sliceStride = 0;
        for (auto [stride, index] : iterators::const_zip(tensor.strides(), indices)) {
            if (index.isNone()) {
                sliceStride = stride;
            }

            baseOffset += stride * index;
        }

        if (static_cast<long>(baseOffset + (values.size() - 1) * sliceStride) > tensor.numel()) {
            throw std::out_of_range("out of range access to tensor storage");
        }

        for (auto [idx, val] : iterators::enumerate(values)) {
            tensor.data_ptr<T>()[baseOffset + idx * sliceStride] = val;
        }

        return tensor;
    }


    /**
     * Creates a 2xN index tensor from two vectors of length N
     * @param from source indices
     * @param to destination indices
     * @return 2xN index tensor
     */
    auto makeIndexTensor(const std::vector<IndexType> &from, const std::vector<IndexType> &to) -> torch::Tensor;

    /**
     * Creates a 2xN index tensor from a vector of edges
     * @param edges edges used to index task nodes
     * @return 2xN index tensor
     */
    auto makeIndexTensor(const EdgeVector &edges) -> torch::Tensor;

    using IndexSpan = std::span<IndexType>;

    /**
     * Returns one row of a COO index tensor
     * @param indexTensor tensor to index
     * @param row row in {0, 1}
     * @return span memory view
     * @throws std::runtime_error if tensor is not contiguous or has the wrong dimensions
     */
    auto getIndexSlice(const torch::Tensor &indexTensor, short row) -> IndexSpan ;

    namespace impl {
        struct MakeEdge {
            constexpr Edge operator()(auto && tuple) const noexcept {
                return std::make_from_tuple<Edge>(std::forward<decltype(tuple)>(tuple));
            }
        };

        using EdgeView = decltype(iterators::zip(std::declval<IndexSpan>(), std::declval<IndexSpan>()) |
                                  std::views::transform(MakeEdge{}));
    }

    /**
     * Creates an iterable view of edges of the edge tensor
     * @param edgeTensor base tensor for the view
     * @return an iterable view that returns Edges
     */
    auto getEdgeView(const torch::Tensor &edgeTensor) -> impl::EdgeView;

    /**
     * Saves a tensor to a file. It can then be loaded from python via torch.load or using loadTensor.
     * @param tensor tensor to save
     * @param fileName destination file
     * @throws std::runtime_error if destination file cannot be opened
     */
    void saveTensor(const torch::Tensor &tensor, const std::filesystem::path &fileName);

    /**
     * Loads a tensor that was saved using saveTensor.
     * @param fileName source file
     * @throws std::runtime_error if source file cannot be opened
     * @return deserialized tensor
     */
    auto loadTensor(const std::filesystem::path &fileName) -> torch::Tensor;

    /**
     * Saves all tensors of a tempo::nn::InputGraph to a directory
     * @param graph graph to save
     * @param fileNameBase path to existing directory
     */
    void saveGraph(const InputGraph &graph, const std::filesystem::path &fileNameBase);

    /**
     * Loads a graph previously saved using saveGraph
     * @param fileNameBase path to directory with tensors
     * @return deserialized graph
     */
    auto loadGraph(const std::filesystem::path &fileNameBase) -> InputGraph;

    /**
     * Checks if two InputGraphs are equivalent. Also prints a detailed report
     * @param graphA first graph
     * @param graphB second graph
     * @return true if the graphs equivalent, i.e. have the same topology and features (permutation invariant), false
     * otherwise
     */
    bool compareGraphs(const InputGraph &graphA, const InputGraph &graphB);

    /**
     * Checks if two InputGraphs are equivalent. Optionally prints a report for the first mismatching attribute.
     * This method is faster than compareGraphs since it stops on the first mismatch
     * @param graphA first graph
     * @param graphB second graph
     * @return true if the graphs equivalent, i.e. have the same topology and features (permutation invariant), false
     * otherwise
     */
    bool graphsEquivalent(const InputGraph &graphA, const InputGraph &graphB, bool verbose = false);
}

#endif //SCHEDCL_TENSOR_UTILS_HPP
