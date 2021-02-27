#include "DirectApproachHR2LQ.h"
#include "BipartiteGraph.h"
#include "StableMarriage.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <climits>
#include <set>
#include <iostream>

DirectApproachHR2LQ::DirectApproachHR2LQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

// runs floyd warshall alg and 
// returns true if negative edge cycle is found
bool DirectApproachHR2LQ::floyd_warshall(std::vector<std::vector<int>>& weight, int num_of_vertices) {
    for (int i = 0; i < num_of_vertices; i++) {
        for (int j = 0; j < num_of_vertices; j++) {
            for (int k = 0; k < num_of_vertices; k++) {
                if (weight[i][j] == INT_MAX || weight[j][k] == INT_MAX) {
                    continue;
                }
                if (weight[i][j] + weight[j][k] < weight[i][k]) {
                    weight[i][k] = weight[i][j] + weight[j][k];
                }
            }
        }
    }

    // checking for cycle
    for (int i = 0; i < num_of_vertices; i++) {
        if (weight[i][i] < 0) {
            return true;
        }
    }
    return false;
}

// runs bellman ford alg and 
// returns true if negative edge cycle is found
bool DirectApproachHR2LQ::bellman_ford(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M, std::map<VertexPtr, int> &index, 
    std::vector<std::vector<int>> &weight, int num_of_vertices) {

    // to maintain distance of source to each vertex
    std::map<VertexPtr, int> distance;

    // VertexBookkeeping for maintaining propose pointer
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // set distances of source to each resident as 0
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;
        distance[r] = 0;
    }

    // set distances of source to each hospital as 0
    for (auto& it : G->get_B_partition()) {
        auto h = it.second;
        distance[h] = 0;
    }

    for (int i = 1; i <= num_of_vertices; i++) {
        //for each edge
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

                // if r->h is an edge
                if (weight[index_of_r][index_of_h] != INT_MAX) {
                    if (distance[r] + weight[index_of_r][index_of_h] < distance[h]) {
                        distance[h] = distance[r] + weight[index_of_r][index_of_h];
                    }
                }

                // if h->r is an edge
                if (weight[index_of_h][index_of_r] != INT_MAX) {
                    if (distance[h] + weight[index_of_h][index_of_r] < distance[r]) {
                        distance[r] = distance[h] + weight[index_of_h][index_of_r];
                    }
                }

                //incrementing propose pointer of r
                bookkeep_data[r].begin += 1;
            }
        }
    }

    // check for negative cycle
    //for each edge
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

            // if r->h is an edge
            if (weight[index_of_r][index_of_h] != INT_MAX) {
                if (distance[r] + weight[index_of_r][index_of_h] < distance[h]) {
                    return true;
                }
            }

            // if h->r is an edge
            if (weight[index_of_h][index_of_r] != INT_MAX) {
                if (distance[h] + weight[index_of_h][index_of_r] < distance[r]) {
                    return true;
                }
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }
    return false;
}

// return true if matching is popular
bool DirectApproachHR2LQ::is_popular(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M) {
    
    // no of vertices will be residents + hospitals
    int num_of_vertices = G->get_A_partition().size() + G->get_B_partition().size();

    // VertexBookkeeping for maintaining propose pointer
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    //to maintain indices of residents and hospitals
    std::map<VertexPtr, int> index;

    // to maintain weights of edges
    std::vector<std::vector<int>> weight(num_of_vertices, std::vector<int>(num_of_vertices, INT_MAX));
    
    int index_count = 0;

    // set indices of residents
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;
        index[r] = index_count;
        weight[index[r]][index[r]] = 0;
        index_count++;
    }
    
    // set indices of hospitals
    for (auto& it : G->get_B_partition()) {
        auto h = it.second;
        index[h] = index_count;
        weight[index[h]][index[h]] = 0;
        index_count++;
    }
    
    // for each edge we will calculate weight
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
        
            bool r_prefers = true;
            bool h_prefers = true;

            // if r has some matchings
            auto M_r = M->find(r);
            if (M_r != M->end()) {
                auto& partners = M_r->second;
                // if edge (r,h) is matched
                // weight will be 0 and direction is from h to r
                if (partners.find(h) != partners.cend()) {
                    weight[index_of_h][index_of_r] = 0;
                    bookkeep_data[r].begin += 1;
                    continue;
                }
                // if (r,h) has unmatched
                else {
                    // if r is fullysubscribed and r prefers least preferred over h
                    if (number_of_partners(M, r) == r->get_upper_quota()) {
                        // r's least preferred partner
                        auto r_least_partner = M->at(r).get_least_preferred().vertex;
                        auto r_least_partner_rank = compute_rank(r_least_partner, r_pref_list);
                        if (r_least_partner_rank < compute_rank(h, r_pref_list)) {
                            r_prefers = false;
                        }
                    }
                }
            }
            else if (r->get_upper_quota() == 0) {
                r_prefers = false;
            }

            if (h->get_upper_quota() == 0) {
                h_prefers = false;
            }
            // if h is fullysubscribed and h prefers least preferred over r
            else if (number_of_partners(M, h) == h->get_upper_quota()) {
                // h's least preferred partner
                auto h_least_partner = M->at(h).get_least_preferred().vertex;
                auto h_least_partner_rank = compute_rank(h_least_partner, h_pref_list);
                if (h_least_partner_rank < compute_rank(r, h_pref_list)) {
                    h_prefers = false;
                }
            }

            // now assign weights
            // if +- or -+, weight = 0
            if ((r_prefers && !h_prefers) || (!r_prefers && h_prefers)) {
                weight[index_of_r][index_of_h] = 0;
            }
            // if ++, weight = -1
            else if (r_prefers && h_prefers) {
                weight[index_of_r][index_of_h] = -1;
            }
            // if --
            else if (!r_prefers && !h_prefers) {
                weight[index_of_r][index_of_h] = 1;
            }
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }
    // call bellman ford algo to check -ve edge cycle
    if (bellman_ford(G, M, index, weight, num_of_vertices)) {
        std::cout << "negative edge cycle found in bellman ford\n";
        return false;
    }
    // call floyd warshall algo to check -ve edge cycle
    if (floyd_warshall(weight, num_of_vertices)) {
        std::cout << "negative edge cycle found in floyd warshall\n";
        return false;
    }
    return true;
}

bool DirectApproachHR2LQ::level_proposing(std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M,
    const BipartiteGraph::ContainerType& proposing_partition) {

    // queue to maintain proposing vertices
    FreeListType free_list;

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // input graph
    std::shared_ptr<BipartiteGraph> G = get_graph();

    // to maintain max level possible of proposing vertices
    unsigned int max_level = 0;

    // set the level of every proposing vertex to 0
    // mark all deficient proposing vetices free (by pushing into the free_list)
    // set max_level to sum of lq's of all proposing vetices
    for (auto& it : proposing_partition) {
        auto v = it.second;
        max_level += v->get_lower_quota();
        // sets level 0
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        // if v is deficient and not present in free list(queue)
        if (number_of_partners(M, v) < v->get_lower_quota()) {
            if (bookkeep_data[v].in_free_list == false) {
                // add it to free list
                bookkeep_data[v].in_free_list = true;
                free_list.push(v);
            }
        }
    }

    // while there is at least one proposing vertex in the free list
    while (not free_list.empty()) {
        // proposing vertex at the front in free list
        auto u = free_list.front();
        //remove u from freelist
        free_list.pop();
        bookkeep_data[u].in_free_list = false;

        auto& u_pref_list = u->get_preference_list();

        // if u hasn't exhausted its preference list
        if (not bookkeep_data[u].is_exhausted()) {
            // highest ranked resident/hospital to whom u has not yet proposed in this level
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            //if v is already fully subscribed
            if (number_of_partners(M, v) == v->get_upper_quota()) {
                // v's least preferred partner
                auto v_partner = M->at(v).get_least_preferred();
                auto possible_partner = Partner(u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                // least preferred partner is compared to present proposing vertex
                // Note : < operator is overloaded
                // if v prefers u over least preferred in its matchings
                if (v_partner < possible_partner) {
                    // remove v_partner from v's matchings and vice versa
                    M->at(v).remove_least_preferred();
                    M->at(v_partner.vertex).remove(v);

                    // add v and u to the matching
                    add_partner(M, u, v, compute_rank(v, u_pref_list), 0);
                    add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                    // add v_partner to free_list if it is not already present
                    if (bookkeep_data[v_partner.vertex].in_free_list == false) {
                        // add it to free list
                        bookkeep_data[v_partner.vertex].in_free_list = true;
                        free_list.push(v_partner.vertex);
                    }
                }
                // if v doesnt prefer u
                else {
                    // add u to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            //if v has vacancies
            else {
                // add v and u to the matching
                add_partner(M, u, v, compute_rank(v, u_pref_list), 0);
                add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
            }
            // if u has vacancies and its level is 0
            if (bookkeep_data[u].level == 0 && number_of_partners(M, u) < u->get_upper_quota()) {
                if (bookkeep_data[u].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            // if u is deficient
            else if (u->get_lower_quota() > 0 && number_of_partners(M, u) < u->get_lower_quota()) {
                if (bookkeep_data[u].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            //increment u's propose pointer
            bookkeep_data[u].begin++;
        }
        // if u is a lq hospital and exhausted its pref list and
        // its level is less than max level
        else if (u->get_lower_quota() > 0 && bookkeep_data[u].level < max_level) {
            // increment h's level
            bookkeep_data[u].level += 1;
            // reset proposal index
            bookkeep_data[u].begin = 0;
            // add it to free list
            bookkeep_data[u].in_free_list = true;
            free_list.push(u);
        }

    }

    // if any proposing vertex is deficient
    // then feasible matching not possible
    for (auto& it : proposing_partition) {
        auto v = it.second;
        if (number_of_partners(M, v) < v->get_lower_quota()) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> DirectApproachHR2LQ::compute_matching() {
    
    std::shared_ptr<BipartiteGraph> G = get_graph();
    
    //find stable matching
    StableMarriage alg(G, false);
    auto M = alg.compute_matching();

    // hospital proposing phase
    if (!level_proposing(M, G->get_B_partition())) {
        std::cout << "Feasible matching not possible\n";
    }

    // resident proposing phase
    if (!level_proposing(M, G->get_A_partition())) {
        std::cout << "Feasible matching not possible\n";
    }
    
    if (!is_popular(G, M)) {
        std::cout << "Not popular\n";
    }

    return M;
}
