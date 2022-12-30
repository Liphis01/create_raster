#ifndef RTREE_H
#define RTREE_H

#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;

typedef vector<pair<double, double>> BR; // Bounding Rectangle

enum NodeType
{
    INTERNAL,
    LEAF,
    DATA
};

/**
 * @brief Struct representing a triangle in 2D space.
 *
 * A triangle is defined by three vertices, each represented as a pair of x and y coordinates.
 */
struct Triangle
{
    pair<double, double> vertex0; ///< (x0, y0)
    pair<double, double> vertex1; ///< (x1, y1)
    pair<double, double> vertex2; ///< (x2, y2)

    /**
     * @brief Constructor for a triangle.
     *
     * Creates a triangle with default values for its vertices.
     */
    Triangle();

    /**
     * @brief Constructor for a triangle.
     *
     * Creates a triangle with the given values for its vertices.
     *
     * @param x0 x coordinate of the first vertex.
     * @param y0 y coordinate of the first vertex.
     * @param x1 x coordinate of the second vertex.
     * @param y1 y coordinate of the second vertex.
     * @param x2 x coordinate of the third vertex.
     * @param y2 y coordinate of the third vertex.
     */
    Triangle(double x0, double y0, double x1, double y1, double x2, double y2);

    /**
     * @brief Check whether a point is inside the triangle.
     *
     * @param point Point to check.
     * @return true if the point is inside the triangle, false otherwise.
     */
    bool pointInTriangle(const pair<double, double> &point) const;

    /**
     * @brief Calculates the area of the triangle.
     *
     * @return Area of the triangle.
     */
    double Area() const;
};

/**
 * @brief Struct representing a node in an R-tree.
 *
 * A node in an R-tree can be either a leaf or an internal node. If it is a leaf, it stores a triangle,
 * otherwise it stores a list of children. Each node also has a bounding box that covers all of its children.
 */
struct RTreeNode
{
    NodeType nodeType;            ///< Whether this node is a leaf or an internal node
    const Triangle dataObject;    ///< Stores the triangle after leaf nodes
    BR bounds;                    ///< Bounding box for this node, covering all its children
    vector<RTreeNode *> children; ///< List of children for this node

    /**
     * @brief Constructor for an R-tree node.
     *
     * @param nodeType Type of the node (leaf or internal).
     */
    RTreeNode(NodeType nodeType);

    /**
     * @brief Constructor for a leaf R-tree node.
     *
     * @param dataObject Triangle to be stored in the leaf node.
     * @param dataBounds Bounding box for the triangle.
     */
    RTreeNode(const Triangle &dataObject, const BR &dataBounds);

    /**
     * @brief Destructor for an R-tree node.
     */
    ~RTreeNode();
};

/**
 * @brief Class representing an R-tree.
 *
 * An R-tree is a tree data structure used for efficient search and insertion of 2D spatial data.
 * It maintains a balanced tree of nodes, each of which has a bounding box covering all of its children.
 * A node can be either a leaf or an internal node. Leaf nodes store triangles, while internal nodes
 * store pointers to other nodes.
 */
class RTree
{
public:
    /**
     * @brief Constructor for an R-tree.
     *
     * Creates an R-tree with default values for the minimum and maximum number of children allowed for each node.
     */
    RTree();

    /**
     * @brief Constructor for an R-tree.
     *
     * Creates an R-tree with the given values for the minimum and maximum number of children allowed for each node.
     *
     * @param MIN_CHILDREN Minimum number of children allowed for each node.
     * @param MAX_CHILDREN Maximum number of children allowed for each node.
     */
    RTree(int MIN_CHILDREN, int MAX_CHILDREN);

    /**
     * @brief Destructor for an R-tree.
     */
    ~RTree();

    /**
     * @brief Inserts a triangle into the R-tree.
     *
     * @param triangle Triangle to be inserted.
     */
    void insert(const Triangle &triangle);

    /**
     * @brief Searches the R-tree for triangles that may contain the given point.
     *
     * @param point Point to search for.
     * @return Vector of triangles that may contain the given point.
     */
    vector<Triangle> search(const pair<double, double> &point) const;

    const int m_m = 2, m_M = 7; ///< Minimum and maximum number of children allowed for each node.
    int m_depth = 1, m_dataNumber = 0; ///< Depth and number of stored data of the tree.

private:
    /**
     * @brief Recursive helper function for searching the R-tree for triangles that may contain the given point.
     *
     * @param point Point to search for.
     * @param node Current node being searched.
     * @return Vector of triangles that may contain the given point.
     */
    vector<Triangle> search(const pair<double, double> &point, const RTreeNode &node) const;
    /**
     * @brief Chooses the leaf node in the R-tree where the given triangle should be inserted.
     *
     * @param triangleBounds Bounding box for the triangle to be inserted.
     * @param path Vector of nodes traversed during the search for the leaf node.
     * @return Pointer to the leaf node where the triangle should be inserted.
     */
    RTreeNode *chooseLeaf(const BR &triangleBounds, vector<RTreeNode *> &path) const;

    /**
     * @brief Splits the given node into two new nodes.
     *
     * @param parentNode Pointer to the parent of the node to be split.
     * @param node Pointer to the node to be split.
     */
    void splitNode(RTreeNode *const parentNode, RTreeNode *node);

    /**
     * @brief Picks the two entries with the largest separation to be the seeds for the split.
     *
     * @param node Pointer to the node whose entries are being split.
     * @param entry1 Pointer to the first seed entry.
     * @param entry2 Pointer to the second seed entry.
     */
    void pickSeeds(const RTreeNode *const node, RTreeNode *&entry1, RTreeNode *&entry2) const;

    /**
     * @brief Picks the next entry for a split, based on the minimum area cost.
     *
     * @param entries Vector of entries being split.
     * @param bounds1 Bounding box for one of the new nodes being created.
     * @param bounds2 Bounding box for the other new node being created.
     * @param areaDiff Difference in area between the two bounding boxes.
     * @return Pointer to the next entry to be added to one of the new nodes.
     */
    RTreeNode *pickNext(const vector<RTreeNode *> &entries, const BR &bounds1, const BR &bounds2, double &areaDiff) const;

    /**
     * @brief Adjusts the tree after an insertion or a split.
     *
     * @param path Vector of nodes traversed during the insertion or split operation.
     * @param newBR Bounding box for the node that was inserted or created during the split.
     */
    void adjustTree(vector<RTreeNode *> &path, const BR &newBR);

    /**
     * @brief Calculates the bounding box for the union of two bounding boxes.
     *
     * @param bounds1 Bounding box for the first region.
     * @param bounds2 Bounding box for the second region.
     * @return Bounding box for the union of the two regions.
     */
    BR calculateBounds(const BR &bounds1, const BR &bounds2) const;

    /**
     * @brief Calculates the area of a bounding box.
     * @param bounds Bounding box for which to calculate the area.
     * @return Area of the bounding box.
    */
    double calculateBoundsArea(const BR &bounds) const;

    /**
        @brief Calculates the cost of adding a bounding box to another bounding box.
        The cost is defined as the increase in the area of the bounding box after adding the new bounding box.
        @param bounds1 Bounding box to which the new bounding box is being added.
        @param boundsOfInsertedNode New bounding box being added.
        @return Cost of adding the new bounding box to the existing bounding box.
        */
    double calculateAreaCost(const BR &bounds1, const BR &boundsOfInsertedNode) const;

    RTreeNode *root_; ///< Pointer to the root node of the tree
};

ostream &operator<<(ostream &os, const Triangle &triangle);
ostream &operator<<(ostream &os, const BR &b);
ostream &operator<<(ostream &os, const RTree &tree);

#endif
