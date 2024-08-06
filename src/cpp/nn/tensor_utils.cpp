/**
 * @author Tim Luchterhand
 * @date 17.07.23.
 */
#include <Iterators.hpp>
#include <cassert>
#include <iostream>
#include <torch/csrc/jit/serialization/pickle.h>
#include <boost/container_hash/hash.hpp>
#include <unordered_set>
#include <ranges>

#include "nn/tensor_utils.hpp"
#include "nn/torch_types.hpp"
#include "util/Matrix.hpp"

namespace tempo::nn::util{
    auto makeIndexTensor(const std::vector<IndexType> &from, const std::vector<IndexType> &to) -> torch::Tensor {
        assert(from.size() == to.size());
        auto ret = torch::empty({2, static_cast<long>(from.size())}, indexTensorOptions());
        for (auto [idx, f, t] : iterators::zip_enumerate(from, to, 0l)) {
            index<IndexType>(ret, {0, idx}) = f;
            index<IndexType>(ret, {1, idx}) = t;
        }

        return ret;
    }

    auto makeIndexTensor(const EdgeVector &edges) -> torch::Tensor {
        auto ret = torch::empty({2, static_cast<long>(edges.size())}, indexTensorOptions());
        for (auto [idx, edge] : iterators::enumerate(edges, 0l)) {
            index<IndexType>(ret, {0, idx}) = edge.first;
            index<IndexType>(ret, {1, idx}) = edge.second;
        }

        return ret;
    }

    auto getIndexSlice(const torch::Tensor &indexTensor, short row) -> IndexSpan {
        assert(row == 0 or row == 1);
        if (indexTensor.sizes().size() != 2 or indexTensor.size(0) != 2) {
            throw std::runtime_error("invalid index tensor");
        }

        if (not indexTensor.is_contiguous() or indexTensor.stride(1) != 1) {
            throw std::runtime_error("non-contiguous tensor cannot be indexed by this function");
        }

        const long offset = row * indexTensor.stride(0);
        return {indexTensor.data_ptr<IndexType>() + offset,
                static_cast<std::span<IndexType>::size_type>(indexTensor.size(1))};
    }

    auto getEdgeView(const torch::Tensor &edgeTensor) -> impl::EdgeView {
        return iterators::zip(getIndexSlice(edgeTensor, 0), getIndexSlice(edgeTensor, 1)) |
               std::views::transform(impl::MakeEdge{});
    }

    void saveTensor(const torch::Tensor &tensor, const std::filesystem::path &fileName) {
        auto bytes = torch::jit::pickle_save(tensor);
        std::ofstream file(fileName, std::ios::binary);
        if (not file.is_open()) {
            throw std::runtime_error("unable to open file " + fileName.string() + " for writing");
        }

        file.write(bytes.data(), static_cast<long>(bytes.size()));
    }

    auto loadTensor(const std::filesystem::path &fileName) -> torch::Tensor {
        std::ifstream file(fileName, std::ios::binary);
        if (not file.is_open()) {
            std::cerr << fileName << std::endl;
            throw std::runtime_error("unable to open file " + fileName.string() + " for reading");
        }

        file >> std::noskipws;
        std::vector<char> buffer{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
        torch::IValue x = torch::jit::pickle_load(buffer);
        return x.toTensor();
    }

    void saveGraph(const InputGraph &graph, const std::filesystem::path &fileNameBase) {
        if (not std::filesystem::is_directory(fileNameBase)) {
            throw std::runtime_error(fileNameBase.string() + " is not a valid directory");
        }

        for (const auto &entry : graph) {
            saveTensor(entry.value(), (fileNameBase / entry.key()).string() + ".pt");
        }
    }

    auto loadGraph(const std::filesystem::path &fileNameBase) -> InputGraph {
        if (not std::filesystem::is_directory(fileNameBase)) {
            throw std::runtime_error(fileNameBase.string() + " is not a valid directory");
        }

        InputGraph ret;
        for (auto key : GraphKeys::AllKeys) {
            auto fileName = (fileNameBase / key).string() + ".pt";
            ret.insert(key, loadTensor(fileName));
        }

        return ret;
    }

    namespace impl {
        using EdgeSet = std::unordered_set<tempo::nn::Edge, boost::hash<tempo::nn::Edge>>;

        template<std::ranges::range R>
        static void printRange(const R &range) {
            std::cout << "[";
            for (const auto &elem : range) {
                std::cout << elem << ", ";
            }

            std::cout << "]" << std::endl;
        }

        template<typename T, typename ...Ts>
        static auto setIntersection(const std::unordered_set<T, Ts...> &a, const std::unordered_set<T, Ts...> &b) {
            std::remove_cvref_t<decltype(a)> intersection;
            std::remove_cvref_t<decltype(a)> diff;
            for (const auto &elemA : a) {
                if (b.contains(elemA)) {
                    intersection.emplace(elemA);
                } else {
                    diff.emplace(elemA);
                }
            }

            for (const auto &elemB : b) {
                if (a.contains(elemB)) {
                    intersection.emplace(elemB);
                } else {
                    diff.emplace(elemB);
                }
            }

            return std::pair(intersection, diff);
        }

        static auto getAdjacency(const torch::Tensor &cooMat) -> Matrix<bool> {
            if (cooMat.numel() == 0) {
                return {};
            }

            if (cooMat.sizes().size() != 2 or cooMat.sizes().front() != 2) {
                throw std::runtime_error("invalid coo tensor");
            }

            auto maxIdx = cooMat.max().item<IndexType>() + 1;
            Matrix<bool> adj(maxIdx, maxIdx);
            for (auto [f, t] : getEdgeView(cooMat)) {
                adj(f, t) = true;
            }

            return adj;
        }

        static bool compareIndexTensors(const torch::Tensor &original, const torch::Tensor &refactored,
                                        std::string_view msg, bool verbose) {
            using namespace tempo::nn;
            if (verbose) {
                std::cout << ">>>>>>>>>>>>>>>>> " << msg << " <<<<<<<<<<<<<<<<<" << std::endl;
            }

            if (original.sizes() != refactored.sizes()) {
                if (verbose) {
                    std::cout << "dimensions not equal: " << original.sizes() << " vs " << refactored.sizes() << std::endl;
                } else {
                    return false;
                }
            }

            const auto adjO = getAdjacency(original);
            const auto adjR = getAdjacency(refactored);

            if (adjR.rawData() == adjO.rawData()) {
                if (verbose) {
                    std::cout << "----- tensors equal -----" << std::endl;
                }

                return true;
            }

            if (verbose) {
                const auto origEdges = util::getEdgeView(original);
                const auto refacEdges = util::getEdgeView(refactored);
                EdgeSet originalEdges(origEdges.begin(), origEdges.end());
                EdgeSet refactoredEdges(refacEdges.begin(), refacEdges.end());
                auto [intersection, diff] = setIntersection(originalEdges, refactoredEdges);
                std::cout << "----- index set intersection -----" << std::endl;
                printRange(intersection);
                std::cout << "----- index set difference -----" << std::endl;
                printRange(diff);
                std::cout << "original size " << originalEdges.size() << std::endl;
                std::cout << "refactored size " << refactoredEdges.size() << std::endl;
                std::cout << "intersection size " << intersection.size() << std::endl;
            }

            return false;
        }

        static bool compareFeatureTensors(const torch::Tensor &original, const torch::Tensor &refactored,
                                          std::string_view msg, bool verbose) {
            using namespace tempo::nn;
            if (verbose) {
                std::cout << ">>>>>>>>>>>>>>>>> " << msg << " <<<<<<<<<<<<<<<<<" << std::endl;
            }

            auto equality = original == refactored;
            bool pass = equality.all().item<bool>();
            if (not verbose) {
                return pass;
            }

            if (not pass) {
                std::cout << "----- original -----" << std::endl;
                std::cout << original << std::endl;
                std::cout << "----- refactored -----" << std::endl;
                std::cout << refactored << std::endl;
            } else {
                std::cout << "tensors match" << std::endl;
            }

            return pass;
        }

        template<std::ranges::random_access_range R>
        static auto getIdx(const R &range, const std::ranges::range_value_t<R> &elem) -> std::optional<long> {
            auto res = std::ranges::find(range, elem);
            if (res == std::ranges::end(range)) {
                return {};
            }

            return static_cast<long>(res - std::ranges::begin(range));
        }

        static bool compareRelationFeatures(const torch::Tensor &original, const torch::Tensor &refactored,
                                            const torch::Tensor &origEdgeIdx, const torch::Tensor &refEdgeIDx,
                                            std::string_view message, bool verbose) {
            using namespace tempo::nn::util;
            if (verbose) {
                std::cout << ">>>>>>>>>>>>>>>>> " << message << "<<<<<<<<<<<<<<<<<" << std::endl;
            }

            if (original.sizes() != refactored.sizes()) {
                if (verbose) {
                    std::cout << "dimensions not equal: " << original.sizes() << " vs " << refactored.sizes() << std::endl;
                }

                return false;
            }

            auto origEdges = getEdgeView(origEdgeIdx);
            auto refEdges = getEdgeView(refEdgeIDx);
            unsigned mismatches = 0;
            for (auto [idx, origEdge] : iterators::enumerate(origEdges, 0l)) {
                auto fOrig = original[idx];
                auto refIdx = getIdx(refEdges, origEdge);
                if (not refIdx.has_value()) {
                    if (verbose) {
                        std::cout << "relation " << origEdge << " does note exist in other relation tensor"
                                  << std::endl;
                        ++mismatches;
                        continue;
                    } else {
                        return false;
                    }
                }

                auto fRef = refactored[*refIdx];
                if (verbose) {
                    std::cout << fOrig << "    ";
                    std::cout << fRef;
                }

                if (not (fOrig == fRef).all().item<bool>()) {
                    if (verbose) {
                        std::cout << " <------------------";
                        ++mismatches;
                    } else {
                        return false;
                    }
                }

                if (verbose) {
                    std::cout << std::endl;
                }
            }

            if (verbose) {
                std::cout << std::endl << mismatches << " mismatches" << std::endl;
            }

            return mismatches == 0;
        }

        static bool compareEdgeResourceRelations(const torch::Tensor &original, const torch::Tensor &refactored,
                                                 const torch::Tensor &origEdgeIdx, const torch::Tensor &refEdgeIDx,
                                                 bool verbose) {
            using namespace tempo::nn::util;
            using namespace std::views;
            if (verbose) {
                std::cout << ">>>>>>>>>>>>>>>>> edge resource relations <<<<<<<<<<<<<<<<<" << std::endl;
            }
            if (original.sizes() != refactored.sizes()) {
                if (verbose) {
                    std::cout << "edge resource relation tensors have different sizes: " << original.sizes() << " vs "
                              << refactored.sizes() << std::endl;
                }

                return false;
            }

            std::vector origEdges(getEdgeView(origEdgeIdx).begin(), getEdgeView(origEdgeIdx).end());
            std::vector refEdges(getEdgeView(refEdgeIDx).begin(), getEdgeView(refEdgeIDx).end());
            std::vector oEdgeResourceRelations(getEdgeView(original).begin(), getEdgeView(original).end());
            std::vector rEdgeResourceRelations(getEdgeView(refactored).begin(), getEdgeView(refactored).end());
            assert(oEdgeResourceRelations.size() == static_cast<std::size_t>(original.size(1)));
            assert(rEdgeResourceRelations.size() == static_cast<std::size_t>(refactored.size(1)));
            unsigned mismatches = 0;
            for (auto [oEdgeIdx, oRes] : oEdgeResourceRelations) {
                auto rEdgeIdx = getIdx(refEdges, origEdges[oEdgeIdx]);
                if (not rEdgeIdx.has_value()) {
                    if (verbose) {
                        std::cout << "relation " << origEdges[oEdgeIdx] << " does note exist in other edge index tensor"
                                  << std::endl;
                        ++mismatches;
                        continue;
                    } else {
                        return false;
                    }
                }

                bool found = false;
                for (auto [rEI, rRI]: rEdgeResourceRelations |
                                      filter([rEdgeIdx](const auto &edge) { return edge.first == *rEdgeIdx; })) {
                    if (rRI == oRes) {
                        found = true;;
                        break;
                    }
                }

                if (not found) {
                    if (verbose) {
                        std::cout << "edge " << origEdges.at(oEdgeIdx) << " relates to resource " << oRes
                                  << " but this relation does not exist in the refactored tensor" << std::endl;
                        ++mismatches;
                    } else {
                        return false;
                    }
                }
            }

            if (verbose) {
                std::cout << std::endl << mismatches << " mismatches" << std::endl;
            }

            return mismatches == 0;
        }

        static auto getEdgesByMask(const torch::Tensor &edges, const torch::Tensor &mask,
                                   tempo::nn::IndexType maskVal) -> EdgeSet {
            using namespace torch::indexing;
            auto edgeSet = edges.index({Slice{None}, mask == maskVal});
            assert(edgeSet.size(0) == 2);
            assert(edgeSet.sizes().size() == 2);
            EdgeSet ret;
            for (auto idx = 0; idx < edgeSet.size(1); ++idx) {
                ret.emplace(edgeSet[0][idx].item<tempo::nn::IndexType>(), edgeSet[1][idx].item<tempo::nn::IndexType>());
            }

            return ret;
        }

        static bool checkEdgePairMask(const torch::Tensor &original, const torch::Tensor &refactored,
                                      const torch::Tensor &origEdgeIdx, const torch::Tensor &refEdgeIDx, bool verbose) {
            using namespace tempo::nn::util;
            if (verbose) {
                std::cout << ">>>>>>>>>>>>>>>>> edge pair mask <<<<<<<<<<<<<<<<<" << std::endl;
            }

            if (original.sizes() != refactored.sizes()) {
                if (verbose) {
                    std::cout << "dimensions not equal: " << original.sizes() << " vs " << refactored.sizes() << std::endl;
                }

                return false;
            }

            std::vector origEdges(getEdgeView(origEdgeIdx).begin(), getEdgeView(origEdgeIdx).end());
            std::vector refEdges(getEdgeView(refEdgeIDx).begin(), getEdgeView(refEdgeIDx).end());
            unsigned mismatches = 0;
            for (auto [oEdgeIdx, edge] : iterators::const_enumerate(origEdges, 0l)) {
                auto rEdgeIdx = getIdx(refEdges, edge);
                if (not rEdgeIdx.has_value()) {
                    if (verbose) {
                        std::cout << "relation " << edge << " does note exist in other edge index tensor" << std::endl;
                        ++mismatches;
                        continue;
                    } else {
                        return false;
                    }
                }

                auto oEdges = getEdgesByMask(origEdgeIdx, original, original[oEdgeIdx].item<tempo::nn::IndexType>());
                auto rEdges = getEdgesByMask(refEdgeIDx, refactored, refactored[*rEdgeIdx].item<tempo::nn::IndexType>());
                if (oEdges != rEdges) {
                    if (verbose) {
                        std::cout << "mismatch for edge group including " << edge << std::endl;
                        printRange(oEdges);
                        std::cout << "vs" << std::endl;
                        printRange(rEdges);
                        ++mismatches;
                    } else {
                        return false;
                    }

                }
            }

            if (verbose) {
                std::cout << std::endl << mismatches << " mismatches" << std::endl;
            }

            return mismatches == 0;
        }
    }

    bool compareGraphs(const InputGraph &graphA, const InputGraph &graphB) {
        using K = tempo::nn::GraphKeys;
        bool result = true;
        const auto &edgesA = graphA.at(K::EdgeIdx);
        const auto &edgesB = graphB.at(K::EdgeIdx);
        result &= impl::compareIndexTensors(edgesA, edgesB, "edge indices", true);
        result &= impl::compareIndexTensors(graphA.at(K::ResourceDependencies), graphB.at(K::ResourceDependencies),
                                            "resource dependencies", true);
        result &= impl::compareEdgeResourceRelations(graphA.at(K::EdgeResourceRelations),
                                                     graphB.at(K::EdgeResourceRelations), edgesA, edgesB, true);
        result &= impl::compareFeatureTensors(graphA.at(K::TaskFeatures), graphB.at(K::TaskFeatures), "task features",
                                              true);
        result &= impl::compareFeatureTensors(graphA.at(K::ResourceFeatures), graphB.at(K::ResourceFeatures),
                                              "resource features", true);
        result &= impl::compareRelationFeatures(graphA.at(K::ResourceConsumptions), graphB.at(K::ResourceConsumptions),
                                                graphA.at(K::ResourceDependencies),
                                                graphB.at(K::ResourceDependencies),
                                                "resource consumptions", true);
        result &= impl::compareRelationFeatures(graphA.at(K::EdgeFeatures), graphB.at(K::EdgeFeatures), edgesA, edgesB,
                                                "edge features", true);
        result &= impl::checkEdgePairMask(graphA.at(K::EdgePairMask), graphB.at(K::EdgePairMask), edgesA, edgesB, true);
        return result;
    }

    bool graphsEquivalent(const InputGraph &graphA, const InputGraph &graphB, bool verbose) {
        using K = tempo::nn::GraphKeys;
        const auto &edgesA = graphA.at(K::EdgeIdx);
        const auto &edgesB = graphB.at(K::EdgeIdx);
        return impl::compareIndexTensors(edgesA, edgesB, "edge indices", verbose) &&
               impl::compareIndexTensors(graphA.at(K::ResourceDependencies), graphB.at(K::ResourceDependencies),
                                         "resource dependencies", verbose) &&
               impl::compareEdgeResourceRelations(graphA.at(K::EdgeResourceRelations),
                                                  graphB.at(K::EdgeResourceRelations), edgesA, edgesB, verbose) &&
               impl::compareFeatureTensors(graphA.at(K::TaskFeatures), graphB.at(K::TaskFeatures), "task features",
                                           verbose) &&
               impl::compareFeatureTensors(graphA.at(K::ResourceFeatures), graphB.at(K::ResourceFeatures),
                                           "resource features", verbose) &&
               impl::compareRelationFeatures(graphA.at(K::ResourceConsumptions), graphB.at(K::ResourceConsumptions),
                                             graphA.at(K::ResourceDependencies), graphB.at(K::ResourceDependencies),
                                             "resource consumptions", verbose) &&
               impl::compareRelationFeatures(graphA.at(K::EdgeFeatures), graphB.at(K::EdgeFeatures), edgesA, edgesB,
                                             "edge features", verbose) &&
               impl::checkEdgePairMask(graphA.at(K::EdgePairMask), graphB.at(K::EdgePairMask), edgesA, edgesB, verbose);
    }
}
