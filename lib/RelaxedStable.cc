#include "RelaxedStable.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <iostream>
#include <set>

RelaxedStable::RelaxedStable(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

//check if output matching is relaxed stable
bool RelaxedStable::is_relaxed_stable(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M) {
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    auto A_partition = G->get_A_partition();

    //for each resident
    //set in_free_list to false for every resident
    // we use in_free_list as in_blocking_pair
    for (auto& it : A_partition) {
        auto v = it.second;
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        bookkeep_data[v].in_free_list = false;
    }
    //for each resident check if it can be in a blocking pair
    for (auto& it : A_partition) {
        //resident
        auto u = it.second;
        auto& u_pref_list = u->get_preference_list();

        // while u hasn't exhausted its preference list
        while (not bookkeep_data[u].is_exhausted()) {
            //hospital to be checked if (r,h) is a blocking pair 
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            //if this hospital is matched to this resident then no need to check for later hospitals
            auto M_u = M->find(u);
            if (M_u != M->end()) {
                auto& partners = M_u->second;
                // v is found in u's matchings
                if (partners.find(v) != partners.cend()) {
                    break;
                }
            }

            //if v is fullysubscribed
            if (number_of_partners(M, v) == v->get_upper_quota()) {
                // v's worst partner
                auto v_worst_partner = M->at(v).get_least_preferred();
                auto v_worst_partner_rank = compute_rank(v_worst_partner.vertex, v_pref_list);
                
                auto u_rank = compute_rank(u, v_pref_list);
                
                //if v prefers its worst partner to u check for next hospital
                if (v_worst_partner_rank < u_rank) {
                    bookkeep_data[u].begin += 1;
                    continue;
                }
            }
            bookkeep_data[u].in_free_list = true;
        }
        // if u is unmatched and in blocking pair
        if (number_of_partners(M, u) == 0 && bookkeep_data[u].in_free_list == true) {
            return false;
        }
    }
    //for each hospital
    auto B_partition = G->get_B_partition();
    for (auto& B1 : B_partition) {
        auto v = B1.second;
        auto v_lower_quota = v->get_lower_quota();

        auto M_v = M->find(v);
        if (M_v != M->end()) {
            //partner list of v
            auto& partners = M_v->second;
            for (const auto& i : partners) {
                auto uc = i.vertex;
                if (bookkeep_data[uc].in_free_list == true) {
                    v_lower_quota = v_lower_quota - 1;
                }
            }
        }
        
        if (v_lower_quota < 0) {
            return false;
        }
    }
    return true;
}

//This function creates a modified graph of input graph where
//for each hospital its upper quota in new graph is set to lower_quota in 
//original graph and each hospital in residents preference list has same rank
//in the new graph.
std::unique_ptr<BipartiteGraph> RelaxedStable::get_modified_graph() {
    //old graph and its partitions
    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto A_partition = G->get_A_partition();
    auto B_partition = G->get_B_partition();

    //partitions of new graph that has to be filled
    BipartiteGraph::ContainerType A, B;

    //for each resident create a new vertex and place it in new graph
    for (auto& A1 : A_partition) {
        A.emplace(A1.first, std::make_shared<Vertex>(A1.first, A1.second->get_lower_quota(), A1.second->get_upper_quota()));
    }
    //for each hospital create a new vertex with upper_quota as lower_quota and place it in new graph
    for (auto& B1 : B_partition) {
        B.emplace(B1.first, std::make_shared<Vertex>(B1.first, B1.second->get_lower_quota(), B1.second->get_lower_quota()));
    }
    //we only need A_Partition in new graph for Classified Popular to work
    //for each resident in new graph create a new preference list
    //All hospitals in the preference list should have same rank for Classified Popular to work
    for (auto& A1 : A_partition) {
        VertexPtr u = A[A1.first];
        PreferenceList& new_pref_list = u->get_preference_list();
        PreferenceList& pref_list = A1.second->get_preference_list();
        for (auto& B1 : pref_list) {
            VertexPtr v = B1.vertex;
            if (new_pref_list.size() == 0) {
                new_pref_list.emplace_back(B[v->get_id()]);
            }
            else {
                //new_pref_list.emplace_back_with_tie(B[v->get_id()]);
            }
        }
    }
    return std::make_unique<BipartiteGraph>(A, B);
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> RelaxedStable::compute_matching() {
    //Input graph and its partitions
    std::shared_ptr<BipartiteGraph> G = get_graph();
    const auto& A_partition = G->get_A_partition();
    const auto& B_partition = G->get_B_partition();

    //For finding minimal feasible matching we need modified graph
    const std::unique_ptr<BipartiteGraph>& G1 = get_modified_graph();
    //ClassifiedPopular alg(G1, true);

    FreeListType free_list;
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    // set the level of every unmatched resident to 1
    // set the level of every matched resident to 0
    //add unmatched residents to free_list
    for (auto& it : A_partition) {
        auto v = it.second;
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        //if resident is unmatched, add it to free list and set level as 1
        if (number_of_partners(M, v) == 0) {
            free_list.push(v);
            bookkeep_data[v].level = 1;
            bookkeep_data[v].in_free_list = true;
        }
    }
    
    // there is at least one resident in the free list
    while (not free_list.empty()) {
        // arbitrary resident in free list
        auto u = free_list.top();
        auto& u_pref_list = u->get_preference_list();
        free_list.pop();

        // if u hasn't exhausted its preference list
        if (not bookkeep_data[u].is_exhausted()) {
            // highest ranked hospital to whom u has not yet proposed
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            
            // if v is under subscribed
            if (number_of_partners(M, v) < v->get_upper_quota()) {
                // add u and v to the matching
                add_partner(M, u, v, (RankType)bookkeep_data[u].begin + 1, 0);
                add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
                bookkeep_data[u].in_free_list = false;
            }
            // if v is fully subscribed
            else {
                //check for level 0 residents in v's partner list
                auto M_v = M->find(v);
                if (M_v != M->end()) {
                    auto& partners = M_v->second;
                    for (const auto& i : partners) {
                        auto uc = i.vertex;
                        //if level 0 resident found
                        if (bookkeep_data[uc].level == 0) {
                            // remove M[uc] from M[v], and M[v] from M[uc]
                            M->at(v).remove(uc);
                            M->at(uc).remove(v);

                            // add u and v to the matching
                            add_partner(M, u, v, (RankType)bookkeep_data[u].begin + 1, 0);
                            add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                            // add uc to free_list
                            free_list.push(uc);
                            bookkeep_data[uc].level = 1;
                            bookkeep_data[uc].in_free_list = true;
                            bookkeep_data[u].in_free_list = false;

                            break;
                        }
                    }
                }
                //if level 0 resident not found
                if(bookkeep_data[u].in_free_list == true) {
                    auto v_worst_partner = M->at(v).get_least_preferred();
                    auto possible_partner = Partner(u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
                    //if v prefers u to its worst partner
                    if (v_worst_partner < possible_partner) {
                        // remove M[v_worst_partner] from M[v], and M[v] from M[v_worst_partner]
                        M->at(v).remove_least_preferred();
                        M->at(v_worst_partner.vertex).remove(v);

                        // add u and v to the matching
                        add_partner(M, u, v, (RankType)bookkeep_data[u].begin + 1, 0);
                        add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                        // add v_worst_partner to free_list
                        free_list.push(v_worst_partner.vertex);
                        bookkeep_data[v_worst_partner.vertex].in_free_list = true;
                        bookkeep_data[u].in_free_list = false;
                    }
                    // if v doesnt prefer, add u again to the list
                    else {
                        free_list.push(u);
                    }
                }
            }
            bookkeep_data[u].begin += 1;
        }
    }
    if (is_relaxed_stable(G, M)) {
        std::cout << "Relaxed_stable\n";
    }
    else {
        std::cout << "Not Relaxed_stable\n";
    }
    return M;
}
