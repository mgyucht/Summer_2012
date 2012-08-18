/*
 * node.cxx
 * --------
 *
 *  node.cxx is the implementation for the functions in the Node class provided
 *  in node.hxx.
 */

#include "node.hxx"

std::vector<double> Node::get_position()
{
  return position;
}

void Node::increment_position(std::vector<double> deltadisp)
{
  for (int i = 0; i < dimensions; i++)
  {
    position[ i ] += deltadisp[ i ];
    delta_pos[ i ] = deltadisp[ i ];
  }
}

std::vector<double> Node::displacement_to(Node node)
{
  std::vector<double> disp(node.dimensions), node_location = node.get_position();
  for (int i = 0; i < node.dimensions; i++)
  {
    disp[i] = node_location[i] - position[i];
  }
  return disp;
}

void Node::set_position(std::vector<double> pos)
{
  for (int i = 0; i < dimensions; i++)
  {
    position[i] = pos[i];
    delta_pos[i] = 0.0;
  }
}

void Node::add_bond(Node node)
{
  spr_connections.append(node);
}

void add_connection(Node node1, Node node2)
{
  node1.add_bond(node2);
  node2.add_bond(node1);
}
