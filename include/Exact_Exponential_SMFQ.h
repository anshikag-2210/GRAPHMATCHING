#ifndef EXACT_EXPONENTIAL_SMFQ_H
#define EXACT_EXPONENTIAL_SMFQ_H

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <vector>
#include <queue>
#include <ostream>

// SEA Popular matching
class Exact_Exponential_SMFQ : public MatchingAlgorithm {
public:
    explicit Exact_Exponential_SMFQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing = false);
    ~Exact_Exponential_SMFQ() override = default;
    void find_indices(std::shared_ptr<BipartiteGraph> G,
        std::map<VertexPtr, unsigned int>& index);
    int gcd(int a, int b);
    int find_costs(std::shared_ptr<BipartiteGraph> G,
        std::map<VertexPtr, unsigned int>& cost);
    void find_distinct_costs(std::shared_ptr<BipartiteGraph> G,
        std::map<VertexPtr, unsigned int>& cost, std::vector<std::vector<int>>& distinct_costs,
        unsigned int &min_cost_possible, unsigned int &max_cost_possible);
    bool find_tuples(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M, std::map<VertexPtr, unsigned int>& index,
        std::vector<std::vector<bool>>& edges, std::map<VertexPtr, unsigned int>& cost,
        std::vector<std::vector<int>>& distinct_costs, std::vector<int>& temp_tuple, int pointer,
        unsigned int required_cost);
    int is_r_perfect(std::vector<int>& degree_of_residents);
    bool find_matching_for_tuple(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> &M, std::map<VertexPtr, unsigned int>& index,
        std::vector<std::vector<bool>>& edges, std::map<VertexPtr, unsigned int>& cost,
        std::vector<int>& temp_tuple);
    bool find_valid_tuples(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M, std::map<VertexPtr, unsigned int>& index,
        std::vector<std::vector<bool>>& edges, std::map<VertexPtr, unsigned int>& cost,
        std::vector<std::vector<std::vector<int>>>& cost_matrix,
        int present_index, int req_cost, std::vector<int>& temp_tuple);
    std::shared_ptr<MatchedPairListType> compute_matching() override;
};

#endif
