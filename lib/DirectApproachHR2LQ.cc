#include "DirectApproachHR2LQ.h"
#include "BipartiteGraph.h"
#include "StableMarriage.h"
#include "Vertex.h"
#include "Partner.h"
#include "Graph.h"
#include "Utils.h"
#include <climits>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>

DirectApproachHR2LQ::DirectApproachHR2LQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

// runs floyd warshall alg and 
// returns true if negative edge cycle is found
bool DirectApproachHR2LQ::floyd_warshall(std::vector<std::vector<int>>& weight, int num_of_vertices,
    std::vector<bool> &lower_quota_vertex) {
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
        if (lower_quota_vertex[i])continue;
        for (int j = 0; j < num_of_vertices; j++) {
            if (lower_quota_vertex[j])continue;
            if (weight[i][j] < 0) {
                std::cout << i << " , " << j << "\n";
                return true;
            }
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
    // for each edge
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

// Given input graph and matching
// calculate weights of all possible edges based on votes of vertices
// call bellman ford to check for negative edge cycle
// if found return false because that matching is not popular
bool DirectApproachHR2LQ::is_popular(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M) {
    
    // no of vertices will be residents + hospitals
    int num_of_vertices = G->get_A_partition().size() + G->get_B_partition().size();

    // VertexBookkeeping for maintaining propose pointer
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    //to maintain indices of residents and hospitals
    std::map<VertexPtr, int> index;

    // to mark lower quota vertices
    std::vector<bool> lower_quota_vertex(num_of_vertices, false);

    // to maintain weights of edges
    std::vector<std::vector<int>> weight(num_of_vertices, std::vector<int>(num_of_vertices, INT_MAX));
    
    int index_count = 0;

    // set indices of residents
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;
        index[r] = index_count;
        weight[index[r]][index[r]] = 0;
        if (r->get_lower_quota() > 0) {
            lower_quota_vertex[index[r]] = true;
        }
        index_count++;
    }
    
    // set indices of hospitals
    for (auto& it : G->get_B_partition()) {
        auto h = it.second;
        index[h] = index_count;
        weight[index[h]][index[h]] = 0;
        if (h->get_lower_quota() > 0) {
            lower_quota_vertex[index[h]] = true;
        }
        index_count++;
    }
    
    // for each edge we will calculate weight based on vote of participants
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
    /*
    // call floyd warshall algo to check -ve edge cycle
    if (floyd_warshall(weight, num_of_vertices, lower_quota_vertex)) {
        std::cout << "negative edge cycle found in floyd warshall\n";
        return false;
    }
    */
    return true;
}

// runs level proposing algorithm
// returns true if feasible matching is possible
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
                //std::cout <<v->get_id()<< " added to free list 1\n";
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

                    //std::cout << v->get_id() << " matched to " << u->get_id() << "\n";

                    // add v_partner to free_list if it is not already present
                    if (bookkeep_data[v_partner.vertex].in_free_list == false) {
                        //std::cout << v_partner.vertex->get_id() << " added to free list 2\n";
                        // add it to free list
                        bookkeep_data[v_partner.vertex].in_free_list = true;
                        free_list.push(v_partner.vertex);
                    }
                }
                // if v doesnt prefer u
                else {
                    // add u to free list
                    //std::cout << u->get_id() << " added to free list 3\n";
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            //if v has vacancies
            else {
                // add v and u to the matching
                add_partner(M, u, v, compute_rank(v, u_pref_list), 0);
                add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
                //std::cout << v->get_id() << " matched to " << u->get_id() << "\n";
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
            if (u->get_lower_quota() > 0 && number_of_partners(M, u) < u->get_lower_quota()) {
                if (bookkeep_data[u].in_free_list == false) {
                    //std::cout << u->get_id() << " added to free list 4\n";
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
            //std::cout << u->get_id() << " added to free list 5 with level "<< bookkeep_data[u].level <<"\n";
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

IdType DirectApproachHR2LQ::get_node_name(IdType id1, IdType id2, IdType id3) {
    return id1 + "_" + id2 + "_" + id3;
}

// To print brandl kavitha graph for SM2LQ instance
void DirectApproachHR2LQ::generate_brandl_kavitha_graph_sm2lq(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M) {

    // to store nodes
    std::vector<NodePtr> nodes;

    // to store the start position of node pointer and number of nodes
    std::map<VertexPtr, std::vector<int>> node_pointers;

    // to store the start position of dummy node pointer and number of nodes for dummies
    std::map<VertexPtr, std::vector<int>> dummynode_pointers;

    // to store the position of unmatched node pointer of any vertex
    std::map<VertexPtr, int> unmatched_node_pointer;

    // to store the position of unmatched dummy node pointer of any vertex
    std::map<VertexPtr, int> unmatched_dummy_node_pointer;

    int present_node_pointer = 0;

    // First find total number of nodes in the graph
    // create those nodes
    
    // for every resident create uq many copies plus uq-lq many dummies
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        // set unmatched node pos to 0
        unmatched_node_pointer[r] = 0;
        // set unmatched dummy node pos to 0
        unmatched_dummy_node_pointer[r] = 0;
        // starting node pointer
        node_pointers[r].push_back(present_node_pointer);
        // number of copies of this vertex
        node_pointers[r].push_back(r->get_upper_quota());
        // create upper quota many copy nodes
        for (int i = 1; i <= r->get_upper_quota(); i++) {
            IdType node_name = get_node_name(r->get_id(), "copy", std::to_string(i));
            NodePtr node(new Node(node_name, present_node_pointer));
            nodes.push_back(node);
            present_node_pointer++;
        }
        // if r has dummies
        if (r->get_upper_quota() - r->get_lower_quota() > 0) {
            // starting dummy node pointer
            dummynode_pointers[r].push_back(present_node_pointer);
            // number of dummies for this vertex
            dummynode_pointers[r].push_back(r->get_upper_quota() - r->get_lower_quota());
            // create uq-lq many dummy nodes
            for (int i = 1; i <= r->get_upper_quota() - r->get_lower_quota(); i++) {
                IdType node_name = get_node_name(r->get_id(), "dummy", std::to_string(i));
                NodePtr node(new Node(node_name, present_node_pointer));
                nodes.push_back(node);
                present_node_pointer++;
            }
        }
    }

    // for every hospital uq many copies plus uq-lq many dummies
    for (auto& it : G->get_B_partition()) {
        //hospital
        auto h = it.second;
        // set unmatched node pos to 0
        unmatched_node_pointer[h] = 0;
        // set unmatched dummy node pos to 0
        unmatched_dummy_node_pointer[h] = 0;
        // starting node pointer
        node_pointers[h].push_back(present_node_pointer);
        // number of copies of this vertex
        node_pointers[h].push_back(h->get_upper_quota());
        // create uq many copy nodes
        for (int i = 1; i <= h->get_upper_quota(); i++) {
            IdType node_name = get_node_name(h->get_id(), "copy", std::to_string(i));
            NodePtr node(new Node(node_name, present_node_pointer));
            nodes.push_back(node);
            present_node_pointer++;
        }
        // if h has dummies
        if (h->get_upper_quota() - h->get_lower_quota() > 0) {
            // starting dummy node pointer
            dummynode_pointers[h].push_back(present_node_pointer);
            // number of dummies for this vertex
            dummynode_pointers[h].push_back(h->get_upper_quota() - h->get_lower_quota());
            // create uq-lq many dummy nodes
            for (int i = 1; i <= h->get_upper_quota() - h->get_lower_quota(); i++) {
                IdType node_name = get_node_name(h->get_id(), "dummy", std::to_string(i));
                NodePtr node(new Node(node_name, present_node_pointer));
                nodes.push_back(node);
                present_node_pointer++;
            }
        }
    }

    // create graph for that many nodes
    Graph cloned_graph = Graph(present_node_pointer);

    // to maintain rank of matched partner for each node
    std::vector<int> matched_node_rank(present_node_pointer, INT_MAX);

    // first consider all matched partners
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        auto M_r = M->find(r);
        // if r has some matchings in M
        if (M_r != M->end()) {
            auto& partners = M_r->second;
            for (const auto& i : partners) {
                // h matched to r in M
                auto h = i.vertex;
                auto& h_pref_list = h->get_preference_list();
                // match some unmatched copy of r to some unmatched copy of h
                // weight will be 0
                int unmatched_pos_of_r = node_pointers[r][0] + unmatched_node_pointer[r];
                int unmatched_pos_of_h = node_pointers[h][0] + unmatched_node_pointer[h];
                // adding edge to graph
                cloned_graph.add_edge(nodes[unmatched_pos_of_r], nodes[unmatched_pos_of_h], 0, true);
                // store rank of matched nodes
                matched_node_rank[unmatched_pos_of_r] = compute_rank(h, r_pref_list);
                matched_node_rank[unmatched_pos_of_h] = compute_rank(r, h_pref_list);
                // increment unmatched pointers
                unmatched_node_pointer[r]++;
                unmatched_node_pointer[h]++;
            }
        }
    }

    // next we create edges for dummies of residents
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        // if r has no dummies
        if (r->get_upper_quota() - r->get_lower_quota() == 0) {
            continue;
        }
        // create a matched edge between unmatched copy of r to unmatched dummy of r
        // create an edge between all copies of r to all dummies of r
        for (int i = node_pointers[r][0]; i < node_pointers[r][0] + node_pointers[r][1]; i++) {
            for (int j = dummynode_pointers[r][0]; j < dummynode_pointers[r][0] + dummynode_pointers[r][1]; j++) {
                int unmatched_pos_of_r = node_pointers[r][0] + unmatched_node_pointer[r];
                int unmatched_dummy_of_r = dummynode_pointers[r][0] + unmatched_dummy_node_pointer[r];
                unmatched_node_pointer[r]++;
                unmatched_dummy_node_pointer[r]++;
                // if this copy and dummy of r both are unmatched
                // create a matched edge between them
                if (i == unmatched_pos_of_r && j == unmatched_dummy_of_r) {
                    cloned_graph.add_edge(nodes[i], nodes[j], 0, true);
                }
                // create a unmatched edge between them
                else {
                    // if copy node of r at pos i is matched to true vertex
                    if (matched_node_rank[i] != INT_MAX) {
                        cloned_graph.add_edge(nodes[i], nodes[j], -1, false);
                    }
                    else {
                        cloned_graph.add_edge(nodes[i], nodes[j], 0, false);
                    }
                }
            }
        }
    }

    // next we create edges for dummies of hospitals
    for (auto& it : G->get_B_partition()) {
        //hospital
        auto h = it.second;
        // if h has no dummies
        if (h->get_upper_quota() - h->get_lower_quota() == 0) {
            continue;
        }
        // create a matched edge between unmatched copy of h to unmatched dummy of h
        // create an edge between all copies of h to all dummies of h
        for (int i = node_pointers[h][0]; i < node_pointers[h][0] + node_pointers[h][1]; i++) {
            for (int j = dummynode_pointers[h][0]; j < dummynode_pointers[h][0] + dummynode_pointers[h][1]; j++) {
                int unmatched_pos_of_h = node_pointers[h][0] + unmatched_node_pointer[h];
                int unmatched_dummy_of_h = dummynode_pointers[h][0] + unmatched_dummy_node_pointer[h];
                unmatched_node_pointer[h]++;
                unmatched_dummy_node_pointer[h]++;
                // if this copy and dummy of h both are unmatched
                // create a matched edge between them
                if (i == unmatched_pos_of_h && j == unmatched_dummy_of_h) {
                    cloned_graph.add_edge(nodes[i], nodes[j], 0, true);
                }
                // create a unmatched edge between them
                else {
                    // if copy node of h at pos i is matched to true vertex
                    if (matched_node_rank[i] != INT_MAX) {
                        cloned_graph.add_edge(nodes[i], nodes[j], -1, false);
                    }
                    else {
                        cloned_graph.add_edge(nodes[i], nodes[j], 0, false);
                    }
                }
            }
        }
    }

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // for each unmatched pair(r,h) we create edges
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto& h_pref_list = h->get_preference_list();

            // if r has some matchings
            auto M_r = M->find(r);
            if (M_r != M->end()) {
                auto& partners = M_r->second;
                // if (r,h) is matched
                // we already created edge in graph
                if (partners.find(h) != partners.cend()) {
                    bookkeep_data[r].begin += 1;
                    continue;
                }
            }
        
            // (r,h) is unmatched, 
            // so create edge between all copies of r to all copies of h
            for (int i = node_pointers[r][0]; i < node_pointers[r][0] + node_pointers[r][1]; i++) {
                for (int j = node_pointers[h][0]; j < node_pointers[h][0] + node_pointers[h][1]; j++) {
                    int wt = 0;
                    // if r's copy node at i is matched to dummy node
                    // r's copy prefers this edge with copy of h over edge with dummy
                    if (matched_node_rank[i] == INT_MAX) {
                        wt += 1;
                    }
                    // if r's copy node at i is matched to true node
                    else {
                        //r's copy prefers its matched node
                        if (matched_node_rank[i] < compute_rank(h, r_pref_list)) {
                            wt -= 1;
                        }
                        else if (matched_node_rank[i] > compute_rank(h, r_pref_list)) {
                            wt += 1;
                        }
                    }
                    // similarly for h
                    // if h's copy node at j is matched to dummy node
                    // h's copy prefers this edge with copy of r over edge with dummys
                    if (matched_node_rank[j] == INT_MAX) {
                        wt += 1;
                    }
                    // if h's copy node at j is matched to true node
                    else {
                        //h's copy prefers its matched node
                        if (matched_node_rank[j] < compute_rank(r, h_pref_list)) {
                            wt -= 1;
                        }
                        else if (matched_node_rank[j] > compute_rank(r, h_pref_list)) {
                            wt += 1;
                        }
                    }

                    // create a unmatched edge
                    cloned_graph.add_edge(nodes[i], nodes[j], wt, false);
                }
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }

    cloned_graph.print_graph();

}

// To print brandl kavitha graph for HR2LQ setting
void DirectApproachHR2LQ::generate_brandl_kavitha_graph_hr2lq(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M) {

    // to store nodes
    std::vector<NodePtr> nodes;

    // to store the start position of node pointer and number of nodes
    std::map<VertexPtr, std::vector<int>> node_pointers;

    // to store the start position of dummy node pointer and number of nodes for dummies
    std::map<VertexPtr, std::vector<int>> dummynode_pointers;

    // to store the position of unmatched node pointer of any vertex
    std::map<VertexPtr, int> unmatched_node_pointer;

    // to store the position of unmatched dummy node pointer of any vertex
    std::map<VertexPtr, int> unmatched_dummy_node_pointer;

    int present_node_pointer = 0;

    // First find total number of nodes in the graph
    // create those nodes
    
    // for every resident create uq many copies plus uq-lq many dummies
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        // set unmatched node pos to 0
        unmatched_node_pointer[r] = 0;
        // set unmatched dummy node pos to 0
        unmatched_dummy_node_pointer[r] = 0;
        // starting node pointer
        node_pointers[r].push_back(present_node_pointer);
        // number of copies of this vertex
        node_pointers[r].push_back(r->get_upper_quota());
        // create upper quota many copy nodes
        for (int i = 1; i <= r->get_upper_quota(); i++) {
            IdType node_name = get_node_name(r->get_id(), "copy", std::to_string(i));
            NodePtr node(new Node(node_name, present_node_pointer));
            nodes.push_back(node);
            present_node_pointer++;
        }
        // if r has dummies
        if (r->get_upper_quota() - r->get_lower_quota() > 0) {
            // starting dummy node pointer
            dummynode_pointers[r].push_back(present_node_pointer);
            // number of dummies for this vertex
            dummynode_pointers[r].push_back(r->get_upper_quota() - r->get_lower_quota());
            // create uq-lq many dummy nodes
            for (int i = 1; i <= r->get_upper_quota() - r->get_lower_quota(); i++) {
                IdType node_name = get_node_name(r->get_id(), "dummy", std::to_string(i));
                NodePtr node(new Node(node_name, present_node_pointer));
                nodes.push_back(node);
                present_node_pointer++;
            }
        }
    }

    // for every hospital uq many copies plus uq-lq many dummies
    for (auto& it : G->get_B_partition()) {
        //hospital
        auto h = it.second;
        // set unmatched node pos to 0
        unmatched_node_pointer[h] = 0;
        // set unmatched dummy node pos to 0
        unmatched_dummy_node_pointer[h] = 0;
        // starting node pointer
        node_pointers[h].push_back(present_node_pointer);
        // number of copies of this vertex
        node_pointers[h].push_back(h->get_upper_quota());
        // create uq many copy nodes
        for (int i = 1; i <= h->get_upper_quota(); i++) {
            IdType node_name = get_node_name(h->get_id(), "copy", std::to_string(i));
            NodePtr node(new Node(node_name, present_node_pointer));
            nodes.push_back(node);
            present_node_pointer++;
        }
        // if h has dummies
        if (h->get_upper_quota() - h->get_lower_quota() > 0) {
            // starting dummy node pointer
            dummynode_pointers[h].push_back(present_node_pointer);
            // number of dummies for this vertex
            dummynode_pointers[h].push_back(h->get_upper_quota() - h->get_lower_quota());
            // create uq-lq many dummy nodes
            for (int i = 1; i <= h->get_upper_quota() - h->get_lower_quota(); i++) {
                IdType node_name = get_node_name(h->get_id(), "dummy", std::to_string(i));
                NodePtr node(new Node(node_name, present_node_pointer));
                nodes.push_back(node);
                present_node_pointer++;
            }
        }
    }

    // create graph for that many nodes
    Graph cloned_graph = Graph(present_node_pointer);

    // to maintain rank of matched partner for each node
    std::vector<int> matched_node_rank(present_node_pointer, INT_MAX);

    // first consider all matched partners
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        auto M_r = M->find(r);
        // if r has some matchings in M
        if (M_r != M->end()) {
            auto& partners = M_r->second;
            for (const auto& i : partners) {
                // h matched to r in M
                auto h = i.vertex;
                auto& h_pref_list = h->get_preference_list();
                // match some unmatched copy of r to some unmatched copy of h
                // weight will be 0
                int unmatched_pos_of_r = node_pointers[r][0] + unmatched_node_pointer[r];
                int unmatched_pos_of_h = node_pointers[h][0] + unmatched_node_pointer[h];
                // adding edge to graph
                cloned_graph.add_edge(nodes[unmatched_pos_of_r], nodes[unmatched_pos_of_h], 0, true);
                // store rank of matched nodes
                matched_node_rank[unmatched_pos_of_r] = compute_rank(h, r_pref_list);
                matched_node_rank[unmatched_pos_of_h] = compute_rank(r, h_pref_list);
                // increment unmatched pointers
                unmatched_node_pointer[r]++;
                unmatched_node_pointer[h]++;
            }
        }
    }

    // next we create edges for dummies of residents
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        // if r has no dummies
        if (r->get_upper_quota() - r->get_lower_quota() == 0) {
            continue;
        }

        // if r is matched to its lower quota many vertices
        // then we create edges between all unmatched copies of r and all its dummies
        if(number_of_partners(M, r) == r->get_lower_quota()){
            int unmatched_pos_of_r = node_pointers[r][0] + unmatched_node_pointer[r];
            // all unmatched copies of r
            while(unmatched_pos_of_r < node_pointers[r][0] + node_pointers[r][1]){
                // all dummies of r
                for (int j = dummynode_pointers[r][0]; j < dummynode_pointers[r][0] + dummynode_pointers[r][1]; j++) {
                    cloned_graph.add_edge(nodes[unmatched_pos_of_r], nodes[j], 0, false);
                }
                unmatched_pos_of_r++;
            }
        }
        // create an edge between all copies of r and all dummies of r
        else{
            // all copies of r
            for (int i = node_pointers[r][0]; i < node_pointers[r][0] + node_pointers[r][1]; i++) {
                // all dummies of r
                for (int j = dummynode_pointers[r][0]; j < dummynode_pointers[r][0] + dummynode_pointers[r][1]; j++) {
                    // if copy node of r at pos i is matched to true vertex
                    if (matched_node_rank[i] != INT_MAX) {
                        cloned_graph.add_edge(nodes[i], nodes[j], -1, false);
                    }
                    else {
                        cloned_graph.add_edge(nodes[i], nodes[j], 0, false);
                    }
                }
            }
        }
    }

    // next we create edges for dummies of hospitals
    for (auto& it : G->get_B_partition()) {
        //hospital
        auto h = it.second;
        // if h has no dummies
        if (h->get_upper_quota() - h->get_lower_quota() == 0) {
            continue;
        }

        // if h is matched to its lower quota many vertices
        // then create an edge between all unmatched copies of h and all its dummies
        if(number_of_partners(M, h) == h->get_lower_quota()){
            int unmatched_pos_of_h = node_pointers[h][0] + unmatched_node_pointer[h];
            // all unmatched copies of h
            while(unmatched_pos_of_h < node_pointers[h][0] + node_pointers[h][1]){
                // all dummies of h
                for (int j = dummynode_pointers[h][0]; j < dummynode_pointers[h][0] + dummynode_pointers[h][1]; j++) {
                    cloned_graph.add_edge(nodes[unmatched_pos_of_h], nodes[j], 0, false);
                }
                unmatched_pos_of_h++;
            }
        }
        // create an edge between all copies of h to all dummies of h
        else{
            // all copies of h
            for (int i = node_pointers[h][0]; i < node_pointers[h][0] + node_pointers[h][1]; i++) {
                // all dummies of h
                for (int j = dummynode_pointers[h][0]; j < dummynode_pointers[h][0] + dummynode_pointers[h][1]; j++) {
                    // if copy node of r at pos i is matched to true vertex
                    if (matched_node_rank[i] != INT_MAX) {
                        cloned_graph.add_edge(nodes[i], nodes[j], -1, false);
                    }
                    else {
                        cloned_graph.add_edge(nodes[i], nodes[j], 0, false);
                    }
                }
            }
        }
    }

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // for each unmatched pair(r,h) we create edges
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto& h_pref_list = h->get_preference_list();

            // if r has some matchings
            auto M_r = M->find(r);
            if (M_r != M->end()) {
                auto& partners = M_r->second;
                // if (r,h) is matched
                // we already created edge in graph
                if (partners.find(h) != partners.cend()) {
                    bookkeep_data[r].begin += 1;
                    continue;
                }
            }
        
            // (r,h) is unmatched, 
            // so create an edge between all copies of r and all copies of h
            for (int i = node_pointers[r][0]; i < node_pointers[r][0] + node_pointers[r][1]; i++) {
                for (int j = node_pointers[h][0]; j < node_pointers[h][0] + node_pointers[h][1]; j++) {
                    int wt = 0;
                    // if r's copy node at i is matched to dummy node
                    // r's copy prefers this edge with copy of h over edge with dummy
                    if (matched_node_rank[i] == INT_MAX) {
                        wt += 1;
                    }
                    // if r's copy node at i is matched to true node
                    else {
                        //r's copy prefers its matched node
                        if (matched_node_rank[i] < compute_rank(h, r_pref_list)) {
                            wt -= 1;
                        }
                        else if (matched_node_rank[i] > compute_rank(h, r_pref_list)) {
                            wt += 1;
                        }
                    }
                    // similarly for h
                    // if h's copy node at j is matched to dummy node
                    // h's copy prefers this edge with copy of r over edge with dummys
                    if (matched_node_rank[j] == INT_MAX) {
                        wt += 1;
                    }
                    // if h's copy node at j is matched to true node
                    else {
                        //h's copy prefers its matched node
                        if (matched_node_rank[j] < compute_rank(r, h_pref_list)) {
                            wt -= 1;
                        }
                        else if (matched_node_rank[j] > compute_rank(r, h_pref_list)) {
                            wt += 1;
                        }
                    }

                    // create a unmatched edge
                    cloned_graph.add_edge(nodes[i], nodes[j], wt, false);
                }
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }

    cloned_graph.print_graph();

}

bool DirectApproachHR2LQ::compare(std::vector<int> A, std::vector<int> B) {
    if (A[0] < B[0])return false;
    else if (A[0] > B[0])return true;
    // if levels are same, consider ranks
    else return (A[1] < B[1]);
}


// converts partition A as SM partition
std::shared_ptr<BipartiteGraph> DirectApproachHR2LQ::reduce_partition_to_SM(BipartiteGraph::ContainerType A,
    BipartiteGraph::ContainerType B, bool is_R_phase) {

    //partitions of new graph
    BipartiteGraph::ContainerType A_new, B_new;

    // to store copy vertices of each vertex
    std::map<VertexPtr, std::vector<VertexPtr>> copy_vertices;

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    int sum_of_lq = 0;
    //find sum of lower quotas of partiton A
    for (auto& A1 : A) {
        auto u = A1.second;
        sum_of_lq += u->get_lower_quota();
    }

    // create replicates of vertices in B partition
    // and add them in new B partition for new graph
    // as their pref lists will be changed
    for (auto& B1 : B) {
        auto u = B1.second;
        B_new.emplace(u->get_id(), std::make_shared<Vertex>(u->get_id(), u->get_cloned_for_id(), u->get_lower_quota(), u->get_upper_quota(), u->is_dummy()));
    }

    // for each vertex in A partition
    for (auto& A1 : A) {
        //vertex
        auto u = A1.second;
        auto& u_pref_list = u->get_preference_list();
        
        // if u is a non-lq vertex,
        // replicate this vertex with same order in its pref list
        if (u->get_lower_quota() == 0) {
            A_new.emplace(u->get_id(), std::make_shared<Vertex>(u->get_id(), u->get_cloned_for_id(), u->get_lower_quota(), u->get_upper_quota(), u->is_dummy()));
            PreferenceList& pref_list = A_new[u->get_id()]->get_preference_list();
            // add vertices in preflist order to preflist of copy of that new vertex
            bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
            while (not bookkeep_data[u].is_exhausted()) {
                //vertex v in u's preflist
                auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
                // add replicate of v in new pref list
                pref_list.emplace_back(B_new[v->get_id()]);
                //incrementing propose pointer of r
                bookkeep_data[u].begin += 1;
            }
            // there will be only level 0 copy for this vertex
            copy_vertices[u].push_back(A_new[u->get_id()]);
        }
        // if u is a lq vertex
        else {
            VertexPtr last_dummy_created;
            // create sum_of_lq+1 many copies for that vertex
            for (int i = 0; i < sum_of_lq+1; i++) {
                IdType node_name = get_node_name(u->get_id(), "level", std::to_string(i));
                A_new.emplace(node_name, std::make_shared<Vertex>(node_name, u->get_cloned_for_id(), 0, 1, false));
                // store copy vertex for future use
                copy_vertices[u].push_back(A_new[node_name]);
                PreferenceList& pref_list = A_new[node_name]->get_preference_list();
                // if it is not level 0 copy
                // add last created dummy to front of its preflist
                if (i != 0) {
                    pref_list.emplace_back(last_dummy_created);
                    // add u's copy in dummy's pref list
                    PreferenceList& dummy_pref_list = last_dummy_created->get_preference_list();
                    dummy_pref_list.emplace_back(A_new[node_name]);
                }
                // add vertices in preflist order to preflist of copy of that vertex
                bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
                while (not bookkeep_data[u].is_exhausted()) {
                    //vertex v in u's preflist
                    auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
                    // add replicate of v in new pref list
                    pref_list.emplace_back(B_new[v->get_id()]);
                    //incrementing propose pointer of r
                    bookkeep_data[u].begin += 1;
                }
                // if it is not last level copy
                // create dummy and add to end of its preflist
                if (i != sum_of_lq) {
                    IdType dummy_node_name = get_node_name(u->get_id(), "dummy", std::to_string(i));
                    B_new.emplace(dummy_node_name, std::make_shared<Vertex>(dummy_node_name, u->get_cloned_for_id(), 0, 1, true));
                    last_dummy_created = B_new[dummy_node_name];
                    pref_list.emplace_back(last_dummy_created);
                    // add u's copy in dummy's pref list
                    PreferenceList& dummy_pref_list = last_dummy_created->get_preference_list();
                    dummy_pref_list.emplace_back(A_new[node_name]);
                }
            }
        }
    }

    // fill preference lists of new B partition
    for (auto& B1 : B) {
        auto u = B1.second;
        auto& u_pref_list = u->get_preference_list();
        PreferenceList& new_pref_list = B_new[u->get_id()]->get_preference_list();
        // for each vertex in u's pref list
        // add their copies in order of levels
        bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
        // to maintain vector of vertices
        std::vector<VertexPtr> vertices;
        int vertex_pointer = 0;
        // to sort the vertices accoring to their levels
        std::vector<std::vector<int>> to_sort;
        while (not bookkeep_data[u].is_exhausted()) {
            //vertex v in u's preflist
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_rank = compute_rank(v, u_pref_list);
            // go through all copies of vertex v
            for (int i = 0; i < copy_vertices[v].size();i++) {
                // add the copy to vector
                vertices.push_back(copy_vertices[v][i]);
                // create a vector containing level, rank, position of this copy
                std::vector<int> temp;
                // level
                temp.push_back(i);
                // rank of v in u's pref_list
                temp.push_back(v_rank);
                // position of this vertex in 'vertices' vector
                temp.push_back(vertex_pointer);

                // add this temp vector to to_sort
                to_sort.push_back(temp);
                vertex_pointer++;
            }
            //incrementing propose pointer of r
            bookkeep_data[u].begin += 1;
        }
        // sort the vertices on basis of level
        sort(to_sort.begin(), to_sort.end(), compare);
        // add those vertices in that order to new pref list
        for (int i = 0; i < to_sort.size(); i++) {
            int pos = to_sort[i][2];
            new_pref_list.emplace_back(vertices[pos]);
        }
    }

    if (is_R_phase) {
        return std::make_shared<BipartiteGraph>(A_new, B_new);
    }
    else {
        return std::make_shared<BipartiteGraph>(B_new, A_new);
    }
}


// make partition A as HR partition
std::shared_ptr<BipartiteGraph> DirectApproachHR2LQ::reduce_partition_to_HR(BipartiteGraph::ContainerType A,
    BipartiteGraph::ContainerType B, bool is_R_phase) {

    //partitions of new graph
    BipartiteGraph::ContainerType A_new, B_new;

    // to store copy vertices of each vertex
    std::map<VertexPtr, std::vector<VertexPtr>> copy_vertices;

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    int sum_of_lq = 0;
    //find sum of lower quotas of partiton A
    for (auto& A1 : A) {
        auto u = A1.second;
        sum_of_lq += u->get_lower_quota();
    }

    // create replicates of vertices in B partition
    // and add them in new B partition for new graph
    // as their pref lists will be changed
    for (auto& B1 : B) {
        auto u = B1.second;
        B_new.emplace(u->get_id(), std::make_shared<Vertex>(u->get_id(), u->get_cloned_for_id(), u->get_lower_quota(), u->get_upper_quota(), u->is_dummy()));
    }

    // for each vertex in A partition
    for (auto& A1 : A) {
        //vertex
        auto u = A1.second;
        auto& u_pref_list = u->get_preference_list();

        // if u is a non-lq vertex,
        // replicate this vertex with same order in its pref list
        if (u->get_lower_quota() == 0) {
            A_new.emplace(u->get_id(), std::make_shared<Vertex>(u->get_id(), u->get_cloned_for_id(), u->get_lower_quota(), u->get_upper_quota(), u->is_dummy()));
            PreferenceList& pref_list = A_new[u->get_id()]->get_preference_list();
            // add vertices in preflist order to preflist of copy of that new vertex
            bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
            while (not bookkeep_data[u].is_exhausted()) {
                //vertex v in u's preflist
                auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
                // add replicate of v in new pref list
                pref_list.emplace_back(B_new[v->get_id()]);
                //incrementing propose pointer of r
                bookkeep_data[u].begin += 1;
            }
            // there will be only level 0 copy for this vertex
            copy_vertices[u].push_back(A_new[u->get_id()]);
        }
        // if u is a lq vertex
        else {
            std::vector<VertexPtr> last_level_dummies;
            // create sum_of_lq+1 many copies for that vertex
            for (int i = 0; i < sum_of_lq + 1; i++) {
                IdType node_name = get_node_name(u->get_id(), "level", std::to_string(i));
                int capacity;
                // level 0 copy; capacity will be upper quota
                if (i == 0) {
                    capacity = u->get_upper_quota();
                }
                // otherwise capacity will be lower quota
                else {
                    capacity = u->get_lower_quota();
                }
                // create the copy vertex
                A_new.emplace(node_name, std::make_shared<Vertex>(node_name, u->get_cloned_for_id(), 0, capacity, false));
                // store copy vertex for future use
                copy_vertices[u].push_back(A_new[node_name]);
                PreferenceList& pref_list = A_new[node_name]->get_preference_list();

                // if it is not level 0 copy
                // add last lower quota many dummies of previous level at the start of pref list
                if (i != 0) {
                    for (int j = last_level_dummies.size() - u->get_lower_quota(); j < last_level_dummies.size(); j++) {
                        pref_list.emplace_back(last_level_dummies[j]);

                        // add u's copy in dummy's pref list
                        PreferenceList& dummy_pref_list = last_level_dummies[j]->get_preference_list();
                        dummy_pref_list.emplace_back(A_new[node_name]);
                    }
                }

                // add vertices in preflist order to preflist of copy of that vertex
                bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
                while (not bookkeep_data[u].is_exhausted()) {
                    //vertex v in u's preflist
                    auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
                    // add replicate of v in new pref list
                    pref_list.emplace_back(B_new[v->get_id()]);
                    //incrementing propose pointer of r
                    bookkeep_data[u].begin += 1;
                }

                // if it is not last level copy
                // we have to create capacity many dummies of that copy and
                // add it to the end of that copies pref list
                if (i != sum_of_lq) {
                    // clear last level dummies to store now created dummies
                    last_level_dummies.clear();
                    for (int j = 0; j < capacity; j++) {
                        IdType dummy_node_name = get_node_name(node_name, "dummy", std::to_string(j));
                        B_new.emplace(dummy_node_name, std::make_shared<Vertex>(dummy_node_name, u->get_cloned_for_id(), 0, 1, true));
                        last_level_dummies.push_back(B_new[dummy_node_name]);

                        // add this dummy to pref list of the vertex
                        pref_list.emplace_back(B_new[dummy_node_name]);

                        // add u's copy in dummy's pref list
                        PreferenceList& dummy_pref_list = B_new[dummy_node_name]->get_preference_list();
                        dummy_pref_list.emplace_back(A_new[node_name]);
                    }
                }
            }
        }
    }

    // fill preference lists of new B partition
    for (auto& B1 : B) {
        auto u = B1.second;
        auto& u_pref_list = u->get_preference_list();
        PreferenceList& new_pref_list = B_new[u->get_id()]->get_preference_list();
        // for each vertex in u's pref list
        // add their copies in order of levels
        bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());
        // to maintain vector of vertices
        std::vector<VertexPtr> vertices;
        int vertex_pointer = 0;
        // to sort the vertices accoring to their levels
        std::vector<std::vector<int>> to_sort;
        while (not bookkeep_data[u].is_exhausted()) {
            //vertex v in u's preflist
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;

            // add all initial dummies to the pref list
            while (v->is_dummy()) {
                new_pref_list.emplace_back(A_new[v->get_id()]);
                //incrementing propose pointer of r
                bookkeep_data[u].begin += 1;
                if (bookkeep_data[u].is_exhausted()) {
                    break;
                }
                v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            }

            if (bookkeep_data[u].is_exhausted()) {
                break;
            }

            // now for all non dummy vertices add them according to their level
            while (!(v->is_dummy())) {
                auto v_rank = compute_rank(v, u_pref_list);
                // go through all copies of vertex v
                for (int i = 0; i < copy_vertices[v].size(); i++) {
                    // add the copy to vector
                    vertices.push_back(copy_vertices[v][i]);
                    // create a vector containing level, rank, position of this copy
                    std::vector<int> temp;
                    // level
                    temp.push_back(i);
                    // rank of v in u's pref_list
                    temp.push_back(v_rank);
                    // position of this vertex in 'vertices' vector
                    temp.push_back(vertex_pointer);

                    // add this temp vector to to_sort
                    to_sort.push_back(temp);
                    vertex_pointer++;
                }
                //incrementing propose pointer of r
                bookkeep_data[u].begin += 1;
                if (bookkeep_data[u].is_exhausted()) {
                    break;
                }
                v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            }
            
            // sort the vertices on basis of level
            sort(to_sort.begin(), to_sort.end(), compare);
            // add those vertices in that order to new pref list
            for (int i = 0; i < to_sort.size(); i++) {
                int pos = to_sort[i][2];
                new_pref_list.emplace_back(vertices[pos]);
            }

            if (bookkeep_data[u].is_exhausted()) {
                break;
            }

            // add all last dummies to the pref list
            while (v->is_dummy()) {
                new_pref_list.emplace_back(A_new[v->get_id()]);
                //incrementing propose pointer of r
                bookkeep_data[u].begin += 1;
                if (bookkeep_data[u].is_exhausted()) {
                    break;
                }
                v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            }
        }
    }

    if (is_R_phase) {
        return std::make_shared<BipartiteGraph>(A_new, B_new);
    }
    else {
        return std::make_shared<BipartiteGraph>(B_new, A_new);
    }
}

std::shared_ptr<BipartiteGraph> DirectApproachHR2LQ::get_reduced_graph() {

    // input graph
    std::shared_ptr<BipartiteGraph> G = get_graph();

    /*
    // hospital phase
    // For finding reduced SM instance graph 
    std::shared_ptr<BipartiteGraph> G1 = reduce_partition_to_SM(G->get_B_partition(), G->get_A_partition(), false);

    // resident phase
    // For finding reduced SM instance graph 
    std::shared_ptr<BipartiteGraph> G2 = reduce_partition_to_SM(G1->get_A_partition(), G1->get_B_partition(), true);
    */

    // hospital phase
    // For finding reduced HR instance graph 
    std::shared_ptr<BipartiteGraph> G1 = reduce_partition_to_HR(G->get_B_partition(), G->get_A_partition(), false);

    // resident phase
    // For finding reduced HR instance graph 
    std::shared_ptr<BipartiteGraph> G2 = reduce_partition_to_HR(G1->get_A_partition(), G1->get_B_partition(), true);

    return G2;
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> DirectApproachHR2LQ::compute_matching() {
    
    std::shared_ptr<BipartiteGraph> G = get_graph();
    
    //below commented code uses level proposing algorithm
    /*
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
    */

    //For finding reduced HR instance graph 
    const std::shared_ptr<BipartiteGraph>& G1 = get_reduced_graph();

    
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    
    // To print reduced instance, uncomment below
    /*
    std::cout << "Reduced graph A partition\n";
    for (auto& it : G1->get_A_partition()) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());
        std::cout << "("<<r->get_lower_quota()<<","<<r->get_upper_quota()<<") "<<r->get_id() << " : ";
        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto& h_pref_list = h->get_preference_list();
            if(h_pref_list.find_index(r) == h_pref_list.size()){
                std::cout<<" not exists ";
            }
            std::cout << h->get_id() <<"("<<compute_rank(h,r_pref_list) << ") , ";
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
        std::cout << " \n";
    }

    std::cout << "\nReduced graph B partition\n";
    for (auto& it : G1->get_B_partition()) {
        //hospital
        auto h = it.second;
        auto& h_pref_list = h->get_preference_list();
        bookkeep_data[h] = VertexBookkeeping(0, h->get_preference_list().size());
        std::cout << "("<<h->get_lower_quota()<<","<<h->get_upper_quota()<<") "<<h->get_id() << " : ";
        // while h hasn't exhausted its preference list
        while (not bookkeep_data[h].is_exhausted()) {
            //resident in h's preflist
            auto r = h_pref_list.at(bookkeep_data[h].begin).vertex;
            auto& r_pref_list = r->get_preference_list();
            if(r_pref_list.find_index(h) == r_pref_list.size()){
                std::cout<<" not exists ";
            }
            std::cout << r->get_id() <<"("<<compute_rank(r,h_pref_list) << ") , ";
            //incrementing propose pointer of r
            bookkeep_data[h].begin += 1;
        }
        std::cout << " \n";
    }
    */

    //finding stable matching for reduced instance
    StableMarriage alg(G1, false);
    auto M = alg.compute_matching();
    
    // convert stable matching of reduced instance to appropriate matching for original graph
    auto A_old = G->get_A_partition();
    auto B_old = G->get_B_partition();
    
    auto N = std::make_shared<MatchingAlgorithm::MatchedPairListType>();
    
    for (auto& it : G1->get_A_partition()) {
        auto r = it.second;
        auto M_r = M->find(r);
        // if resident is a dummy
        if (r->is_dummy() == true) {
            continue;
        }
        
        // r in original graph 
        auto r1 = A_old[r->get_cloned_for_id()];
        auto& r1_pref_list = r1->get_preference_list();
     
        if (M_r != M->end()) {
            auto& partners = M_r->second;
            for (const auto& i : partners) {
                auto h = i.vertex;
                // if hospital is a dummy
                if (h->is_dummy() == true) {
                    continue;
                }
                
                // h in original graph 
                auto h1 = B_old[h->get_cloned_for_id()];
                auto& h1_pref_list = h1->get_preference_list();
                
                // add r and h to the matching
                add_partner(N, r1, h1, compute_rank(h, r1_pref_list), 0);
                add_partner(N, h1, r1, compute_rank(r, h1_pref_list), 0);
            }
        }
      
    }
    
    //To print brandl kavitha graph
    // generate_brandl_kavitha_graph_hr2lq(G, N);

    return N;
}
