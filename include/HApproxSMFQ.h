#ifndef H_APPROX_SMFQ_H
#define H_APPROX_SMFQ_H

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <vector>
#include <queue>
#include <ostream>

// SEA Popular matching
class HApproxSMFQ : public MatchingAlgorithm {
public:
    explicit HApproxSMFQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing = false);
    ~HApproxSMFQ() override = default;
    bool find_costs(std::shared_ptr<BipartiteGraph> G,
        std::map<VertexPtr, unsigned int>& cost);
    bool is_r_perfect(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M);
    std::shared_ptr<MatchedPairListType> compute_matching() override;
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> compute_matching1(unsigned int& C);
};

#endif
