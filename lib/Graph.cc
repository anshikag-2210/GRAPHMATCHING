#include "Graph.h"

#include <queue>
#include <iostream>

Graph::Graph()
{}

Graph::Graph(NodeType num_nodes)
	: num_nodes_(num_nodes), edges_(0), adj_list_(num_nodes)
{}

Graph::~Graph()
{}

void Graph::add_edge(NodePtr U, NodePtr V, int weight /* = 0 */, bool matched /* = false */) {
	if (U->get_id() != V->get_id()) {
		// Add (u,v) edge with given capacity and rank
		edges_.push_back(Edge(U, V, weight, matched));

		// Add pointer of this edge in adjacency list of both u and v
		adj_list_[U->get_id()].push_back(edges_.size() - 1);
		adj_list_[V->get_id()].push_back(edges_.size() - 1);
	}
}

void Graph::print_graph() {
	for (int i = 0; i < edges_.size(); i++) {
		std::cout << edges_[i].get_U()->get_name()<<",";
		std::cout << edges_[i].get_V()->get_name() << ",";
		std::cout << edges_[i].get_weight()<< ",";
		std::cout << edges_[i].is_matched()<<"\n";
		//if(edges_[i].is_matched())std::cout << "matched\n";
		//else std::cout << "unmatched\n";
	}
}
