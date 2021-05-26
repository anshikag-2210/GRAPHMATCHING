#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <vector>

#include "Edge.h"
#include "TDefs.h"

// undirected weighted graph.  
class Graph {
public:
	typedef std::vector<Edge> EdgeList;

private:
	NodeType num_nodes_;                            // number of nodes in the graph
	EdgeList edges_;                                // the list of edges
	std::vector<std::vector<NodeType> > adj_list_;  // the adjacency list of the graph 
	
public:
	Graph();
	Graph(NodeType num_nodes);
	~Graph();

	// Adds an edge to the graph
	void add_edge(NodePtr U, NodePtr V, int weight = 0, bool matched = false);

	// print the graph
	void print_graph();
};

#endif