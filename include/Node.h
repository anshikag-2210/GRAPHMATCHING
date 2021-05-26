#ifndef NODE_H
#define NODE_H

#include "TDefs.h"

// Node in a Graph.
class Node {
private:
	IdType node_name_;            // name of the node (for instance, C_a_b)
	NodeType node_id_;            // id of this node
	
public:
	Node();
	Node(IdType node_name, NodeType node_id);
	~Node();
	const IdType& get_name() const;
	const NodeType& get_id() const;
};

#endif