#ifndef EDGE_H
#define EDGE_H

#include "Node.h"
#include "TDefs.h"

// Edge between Nodes in a Graph.
class Edge{
private:
	NodePtr U_, V_;      // nodes which contain the undirected edge (U_- V_)
	int weight_;  // capacity of the edge
	bool matched_; // if this edge is matched or not
public:
	Edge();
	Edge(NodePtr U, NodePtr V, int weight = 0, bool matched = false);
	~Edge();
	const NodePtr& get_U() const;
	const NodePtr& get_V() const;
	const int& get_weight() const;
	void set_weight(int weight);
	bool is_matched();
	void set_matched(bool matched);
};

#endif