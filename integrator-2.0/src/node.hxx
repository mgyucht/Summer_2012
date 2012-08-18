/*
 * node.hxx
 * --------
 *
 *  Provides the node data structure. A node is the smallest component of a
 *  network. It has a position and a radius, and it is connected to a variable
 *  number of other nodes in the network.
 */

#include <vector>

static int current_node = 0;

class Node {

  public:

    /*
     * The node constructor. It only initializes the position of the node.
     */
    Node (std::vector<double> ppos)
    {
      dimensions = ppos.size();
      set_position(ppos);
      node_identifier = current_node;
      current_node++;
    }

    /*
     * get_position returns the current position of the node.
     */
    std::vector<double> get_position();

    /*
     * increment_position takes as its argument the new position of the node
     * relative to its present position.
     */
    void increment_position(std::vector<double> /* deltadisp */);

    /*
     * displacement_to returns the displacement between the present node and the node
     * passed as the argument to the function as a vector.
     */
    std::vector<double> displacement_to(Node /* node */);

    /*
     * add_bond adds a single bond between the current node and the node passed
     * as an argument.
     */
    void add_bond(Node /* node */);

  private:
    unsigned int node_identifier;
    static int dimensions;
    static double radius;
    std::vector<double> position;
    std::vector<double> delta_pos;
    std::vector<Node> spr_connections;
    std::vector<Node> spr_motors;

    /*
     * set_position takes as its argument a vector of doubles corresponding to
     * the location of the node and sets the position of the node to the
     * corresponding values. This function should only be called once to
     * initialize the node.
     */
    void set_position(std::vector<double> /* pos */);

};

/*
 * add_connection creates a new bond between two nodes.
 */
void add_connection(Node /* node1 */, Node /* node2 */);
void remove_connection(Node /* node1 */, Node /* node2 */);
void add_motor(Node /* node1 */, Node /* node2 */);
void remove_motor(Node /* node1 */, Node /* node2 */);
