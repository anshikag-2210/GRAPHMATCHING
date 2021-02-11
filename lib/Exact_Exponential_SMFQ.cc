#include "Exact_Exponential_SMFQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include "StableMarriage.h"
#include "Statistics.h"
#include <set>
#include <climits>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>

Exact_Exponential_SMFQ::Exact_Exponential_SMFQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

// to find indices of residents and hospitals
// indices start from 0 for both residents and hospitals
void Exact_Exponential_SMFQ::find_indices(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& index) {

    int index_count = 0;

    // find indices of residents
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        index[u] = index_count;
        index_count++;
    }

    index_count = 0;

    // find indices of hospitals
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        index[v] = index_count;
        index_count++;
    }
}

// GCD of 'a' and 'b' 
int Exact_Exponential_SMFQ::gcd(int a, int b)
{
    if (b == 0)return a;
    return gcd(b, a % b);
}

// to find costs of hospitals
int Exact_Exponential_SMFQ::find_costs(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int> &cost) {

    // if any resident has empty pref list return -1
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        // if u's pref list is empty
        if (u->get_preference_list().size() == 0) {
            std::cout << "Resident " << u->get_id() << " has no preference list\n";
            return -1;
        }
    }

    unsigned long int lcm = 0;

    // finding lcm's of upper quota's of all hospitals
    for (auto& B1 : G->get_B_partition()) {
        auto v = B1.second;
        // if v's pref list is empty
        if (v->get_preference_list().size() == 0) {
            std::cout << "Hospital " << v->get_id() << " has no preference list\n";
            return -1;
        }
        auto v_upper_quota = v->get_upper_quota();
        //for first hospital
        if (lcm == 0) {
            lcm = v_upper_quota;
        }
        else {
            lcm = (((v_upper_quota * lcm)) / (gcd(v_upper_quota, lcm)));
        }
    }

    // set the cost of every hospital
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        cost[v] = (lcm / (v->get_upper_quota()));
    }
    return 1;
}

void Exact_Exponential_SMFQ::find_distinct_costs(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost,
    std::vector<std::vector<int>>& distinct_costs) {

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // for each resident find all distinct costs in its pref list
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        std::vector<int> distinct_costs_of_r;

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            int cost_of_h = cost[h];

            // if cost(h) is not already present in vector 
            bool found = false;
            for (int p = 0; p < distinct_costs_of_r.size(); p++) {
                if (distinct_costs_of_r[p] == cost_of_h) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                distinct_costs_of_r.push_back(cost_of_h);
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
        //pushing vector to main vector
        distinct_costs.push_back(distinct_costs_of_r);
    }
}

// returns true if edges is r_perfect
int Exact_Exponential_SMFQ::is_r_perfect(std::vector<std::vector<bool>>& edges) {

    for (int i = 0; i < edges.size(); i++) {
        bool found_true = false;
        for (int j = 0; j < edges[i].size(); j++) {
            if (edges[i][j] == true) {
                found_true = true;
                break;
            }
        }
        if (!found_true)return 0;
    }
    return 1;
}

void Exact_Exponential_SMFQ::find_matching_for_tuple(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> &M, std::map<VertexPtr, unsigned int>& index,
    std::vector<std::vector<bool>>& edges, std::map<VertexPtr, unsigned int>& cost, 
    std::vector<int>& temp_tuple, unsigned int &min_cost) {

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // for each resident mark all hospitals in its pref list with respective cost in tuple
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto index_of_r = index[r];
        auto& r_pref_list = r->get_preference_list();

        auto cost_to_be_set = temp_tuple[index_of_r];

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto index_of_h = index[h];

            // if cost[h] is equal to cost_to_be_set
            // we should include edges(r,h) in graph
            if (cost[h] == cost_to_be_set) {
                edges[index_of_r][index_of_h] = true;
            }
            else {
                edges[index_of_r][index_of_h] = false;
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }

    bool change = true;

    // while graph is R-perfect and change is true
    while (is_r_perfect(edges) && change) {

        change = false;

        // for each resident
        for (auto& it : G->get_A_partition()) {
            //resident
            auto r = it.second;
            auto index_of_r = index[r];
            auto& r_pref_list = r->get_preference_list();

            // to find top preferred hospital of r in new pruned graph
            // initially taken as first hospital
            auto top_preffered_hospital_of_r = r_pref_list.at(0).vertex;

            bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

            // while r hasn't exhausted its preference list
            while (not bookkeep_data[r].is_exhausted()) {
                //hospital in r's preflist
                auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
                auto& h_pref_list = h->get_preference_list();
                auto index_of_h = index[h];

                // if hospital h has edge with r in new graph
                if (edges[index_of_r][index_of_h] == true) {
                    top_preffered_hospital_of_r = h;
                    break;
                }
                // h is more preferred by r than top preferred in original graph
                else {
                    // we have to delete all edges (r1, h)
                    // where h prefers r over r1
                    bookkeep_data[h] = VertexBookkeeping(0, h->get_preference_list().size());

                    bool start_deleting = false;
                    // while h hasn't exhausted its preference list
                    while (not bookkeep_data[h].is_exhausted()) {
                        //resident in h's preflist
                        auto r1 = h_pref_list.at(bookkeep_data[h].begin).vertex;
                        auto index_of_r1 = index[r1];

                        // r1 is less preferred than r by h
                        // and edge (r1, h) is present in new graph 
                        if (start_deleting && edges[index_of_r1][index_of_h] == true) {
                            // delete the  edge (r1, h)
                            edges[index_of_r1][index_of_h] = false;
                            change = true;
                        }

                        // if r is occurred in h's preference list
                        // then all next residents are less preferred than r
                        if (index_of_r1 == index_of_r) {
                            start_deleting = true;
                        }

                        //incrementing propose pointer of h
                        bookkeep_data[h].begin += 1;
                    }
                }

                //incrementing propose pointer of r
                bookkeep_data[r].begin += 1;
            }

            // if r exhausted its list and 
            // there is no edge for r in new graph
            auto index_of_top_pref_h = index[top_preffered_hospital_of_r];
            if (edges[index_of_r][index_of_top_pref_h] == false) {
                // this cost tuple wont give r-perfect stable matching
                break;
            }

        } // end of for loop of residents
    }//end of while loop

    // if r-perfect matching is possible
    if (is_r_perfect(edges)) {
        // to maintain matching for this tuple
        auto M1 = std::make_shared<MatchingAlgorithm::MatchedPairListType>();
        auto temp_cost = 0;

        // we will match each resident with
        // first hospital in its preference list
        // which it has edge with in pruned graph
        for (auto& it : G->get_A_partition()) {
            //resident
            auto r = it.second;
            auto index_of_r = index[r];
            auto& r_pref_list = r->get_preference_list();

            bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

            // while r hasn't exhausted its preference list
            while (not bookkeep_data[r].is_exhausted()) {
                //hospital in r's preflist
                auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
                auto index_of_h = index[h];
                auto& h_pref_list = h->get_preference_list();

                // if edge (r,h) is in pruned graph
                if (edges[index_of_r][index_of_h] == true) {
                    // add r and h to the matching M1
                    add_partner(M1, h, r, compute_rank(r, h_pref_list), 0);
                    add_partner(M1, r, h, compute_rank(h, r_pref_list), 0);
                    // add h's cost to temp_cost
                    temp_cost += cost[h];
                    break;
                }

                //incrementing propose pointer of r
                bookkeep_data[r].begin += 1;
            }
        }

        // if M1 is the new min_cost matching
        if (temp_cost < min_cost) {
            min_cost = temp_cost;
            M = M1;
        }
    }
}


// to generate all tuples
void Exact_Exponential_SMFQ::find_tuples(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M, std::map<VertexPtr, unsigned int>& index,
    std::vector<std::vector<bool>>& edges, std::map<VertexPtr, unsigned int>& cost, 
    std::vector<std::vector<int>>& distinct_costs, std::vector<int> &temp_tuple, int pointer) {

    unsigned min_cost = UINT_MAX; 

    //if pointer equal size of distinct costs
    if (pointer >= distinct_costs.size()) {
        find_matching_for_tuple(G, M, index, edges, cost, temp_tuple, min_cost);
        return;
    }
    else {
        // add each cost of that resident one by one
        // and recurse
        for (int i = 0; i < distinct_costs[pointer].size(); i++) {
            temp_tuple.push_back(distinct_costs[pointer][i]);
            find_tuples(G, M, index, edges, cost, distinct_costs, temp_tuple, pointer+1);
            temp_tuple.pop_back();
        }
    }
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> Exact_Exponential_SMFQ::compute_matching() {
    
    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();
    unsigned int R, H;
    R = G->get_A_partition().size();
    H = G->get_B_partition().size();

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    //to maintain indices of residents and hospitals
    std::map<VertexPtr, unsigned int> index;

    //to maintain costs of hospitals
    std::map<VertexPtr, unsigned int> cost;

    //to maintain distinct costs in each residents pref list
    std::vector<std::vector<int>> distinct_costs;

    // to maintain edges of pruned graph for every tuple
    std::vector<std::vector<bool>> edges(R, std::vector<bool>(H,false));

    //find cost of hospitals
    //if hospital or resident has no pref list, find_costs method returns -1
    if (!find_costs(G, cost)) {
        return M;
    }
    
    // find indices of vertices
    find_indices(G, index);

    unsigned long long int count = 1;
    
    // find distinct costs for each resident
    find_distinct_costs(G, cost, distinct_costs);
    
    // find all possible tuples of costs
    std::vector<int> temp_tuple;
    find_tuples(G, M, index, edges, cost, distinct_costs, temp_tuple, 0);

    return M;

}
