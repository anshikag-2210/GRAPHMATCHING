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
    bool level_proposing(std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M,
        const BipartiteGraph::ContainerType& proposing_partition);
    std::shared_ptr<MatchedPairListType> compute_matching() override;
};

#endif
