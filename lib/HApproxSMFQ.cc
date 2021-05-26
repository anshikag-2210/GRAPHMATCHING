#include "HApproxSMFQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include "StableMarriage.h"
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <climits>
#include <iostream>

HApproxSMFQ::HApproxSMFQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

// to find costs of hospitals
// and to print input parameters
bool HApproxSMFQ::find_costs(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost) {

    // find if any resident has empty pref list
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        // if u's pref list is empty
        if (u->get_preference_list().size() == 0) {
            std::cout << "Resident " << u->get_id() << " has no preference list\n";
            return false;
        }
    }

    std::vector<int> capacities;

    // finding maximum of upper quota's of all hospitals
    for (auto& B1 : G->get_B_partition()) {
        auto v = B1.second;
        // if v's pref list is empty
        if (v->get_preference_list().size() == 0) {
            std::cout << "Hospital " << v->get_id() << " has no preference list\n";
            return false;
        }
        
        int v_upper_quota = v->get_upper_quota();
        if (std::find(capacities.begin(), capacities.end(), v_upper_quota) == capacities.end()) {
            capacities.push_back(v_upper_quota);
        }
    }

    // sort the distinct_costs_of_r vector
    sort(capacities.begin(), capacities.end(), std::greater<int>());

    // set the cost of every hospital
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        int v_upper_quota = v->get_upper_quota();
        cost[v] = std::find(capacities.begin(), capacities.end(), v_upper_quota) - capacities.begin() + 1;
    }
   
    return true;
}

bool HApproxSMFQ::is_r_perfect(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M) {
    // check if each resident is matched
    for (auto& it : G->get_A_partition()) {
        auto v = it.second;
        if (number_of_partners(M, v) != v->get_upper_quota()) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> HApproxSMFQ::compute_matching() {

    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    //to maintain costs of hospitals
    std::map<VertexPtr, unsigned int> cost;

    //to save upper quotas of hospitals
    std::map<VertexPtr, unsigned> upper_quotas;

    //find cost of hospitals
    //if hospital or resident has no pref list, find_costs method returns -1
    if (!find_costs(G, cost)) {
        return M;
    }

    // find max_cost of any hospital
    // save upper quota of hospitals
    unsigned int max_cost = 0;
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        if (cost[v] > max_cost)max_cost = cost[v];
        upper_quotas[v] = v->get_upper_quota();
    }

    int a = 0;
    int b = G->get_A_partition().size() * max_cost;

    while (true) {
        int C = floor((a + b)/2);
        // for every hospital h
        // set upper quota as C/cost(h)
        for (auto& it : G->get_B_partition()) {
            auto v = it.second;
            v->set_upper_quota(floor(C/cost[v]));
        }

        // find stable matching for abov instance
        StableMarriage alg(G, true);
        auto Ms = alg.compute_matching();

        // if Ms is R-perfect
        if (is_r_perfect(G, Ms)) {
            // for every hospital h
            // set upper quota as (C-1)/cost(h)
            for (auto& it : G->get_B_partition()) {
                auto v = it.second;
                v->set_upper_quota(floor((C - 1) / cost[v]));
            }

            // find stable matching for above instance
            StableMarriage alg1(G, true);
            auto Ms1 = alg.compute_matching();
            if (is_r_perfect(G, Ms1)) {
                b = C - 1;
            }
            else {
                M = Ms;
                break;
            }
        }
        else {
            a = C + 1;
        }
    }

    // reset upper quotas of all hospitals
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        v->set_upper_quota(upper_quotas[v]);
    }

    return M;
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> HApproxSMFQ::compute_matching1(unsigned int &lb) {

    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    //to maintain costs of hospitals
    std::map<VertexPtr, unsigned int> cost;

    //to save upper quotas of hospitals
    std::map<VertexPtr, unsigned> upper_quotas;

    //find cost of hospitals
    //if hospital or resident has no pref list, find_costs method returns -1
    if (!find_costs(G, cost)) {
        return M;
    }

    // find max_cost of any hospital
    // save upper quota of hospitals
    unsigned int max_cost = 0;
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        if (cost[v] > max_cost)max_cost = cost[v];
        upper_quotas[v] = v->get_upper_quota();
    }

    int a = 0;
    int b = G->get_A_partition().size() * max_cost;

    while (true) {
        int C = floor((a + b) / 2);
        // for every hospital h
        // set upper quota as C/cost(h)
        for (auto& it : G->get_B_partition()) {
            auto v = it.second;
            if (cost[v] != 0) {
                if (floor(C / cost[v]) > v->get_preference_list().size()) {
                    v->set_upper_quota(v->get_preference_list().size());
                }
                else {
                    v->set_upper_quota(floor(C / cost[v]));
                }
            }
            else {
                v->set_upper_quota(v->get_preference_list().size());
            }
        }

        // find stable matching for abov instance
        StableMarriage alg(G, true);
        auto Ms = alg.compute_matching();

        // if Ms is R-perfect
        if (is_r_perfect(G, Ms)) {
            // for every hospital h
            // set upper quota as (C-1)/cost(h)
            for (auto& it : G->get_B_partition()) {
                auto v = it.second;
                if (cost[v] != 0) {
                    if (floor((C-1) / cost[v]) > v->get_preference_list().size()) {
                        v->set_upper_quota(v->get_preference_list().size());
                    }
                    else {
                        v->set_upper_quota(floor((C-1) / cost[v]));
                    }
                }
                else {
                    v->set_upper_quota(v->get_preference_list().size());
                }
            }

            // find stable matching for above instance
            StableMarriage alg1(G, true);
            auto Ms1 = alg.compute_matching();
            if (is_r_perfect(G, Ms1)) {
                b = C - 1;
            }
            else {
                M = Ms;
                lb = C;
                break;
            }
        }
        else {
            a = C + 1;
        }
    }

    // reset upper quotas of all hospitals
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        v->set_upper_quota(upper_quotas[v]);
    }

    return M;
}
