#ifndef LP_APPROX_SMFQ_H
#define LP_APPROX_SMFQ_H

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <vector>
#include <queue>
#include <ostream>

// SEA Popular matching
class LpApproxSMFQ : public MatchingAlgorithm {
public:
    explicit LpApproxSMFQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing = false);
    ~LpApproxSMFQ() override = default;
    int gcd(int a, int b);
    int find_costs(std::shared_ptr<BipartiteGraph> G,
        std::map<VertexPtr, unsigned int>& cost,
        std::vector<std::vector<int>> &additional_output);
    std::shared_ptr<MatchedPairListType> compute_matching() override;
    void print_additional_output(std::shared_ptr<BipartiteGraph> G
        , std::vector<std::vector<int>>& additional_output
        , std::vector<std::string>& additional_output_names);
};

#endif
