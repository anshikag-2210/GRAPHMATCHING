#ifndef DIRECT_APPROACH_HR2LQ_H
#define DIRECT_APPROACH_HR2LQ_H

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <queue>
#include <ostream>

// HR2LQ Direct Approach
class DirectApproachHR2LQ : public MatchingAlgorithm {
private:
    // typedef queue for hospitals
    typedef std::queue<VertexPtr> FreeListType;

public:
    explicit DirectApproachHR2LQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing = false);
    ~DirectApproachHR2LQ() override = default;
    bool floyd_warshall(std::vector<std::vector<int>>& weight, int num_of_vertices, 
        std::vector<bool>& lower_quota_vertex);
    bool bellman_ford(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M, std::map<VertexPtr, int>& index,
        std::vector<std::vector<int>>& weight, int num_of_vertices);
    bool level_proposing(std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M,
        const BipartiteGraph::ContainerType& proposing_partition);
    bool is_popular(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M);
    std::shared_ptr<MatchedPairListType> compute_matching() override;
};

#endif
