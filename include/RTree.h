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

struct Triangle
{
    pair<double, double> vertex0; // (x0, y0)
    pair<double, double> vertex1; // (x1, y1)
    pair<double, double> vertex2; // (x2, y2)

    // Constructor for a triangle
    Triangle();
    Triangle(double, double, double, double, double, double);
    bool PointInTriangle(pair<double, double> &point) const;
    double Area() const;
};

struct RTreeNode
{
    NodeType nodeType;            // Whether this node is a leaf or an internal node
    const Triangle dataObject;          // Stores the triangle after leaf nodes
    BR bounds;                    // Bounding box for this node, covering all its children
    vector<RTreeNode *> children; // List of children for this node

    // Constructor for an R-tree node
    RTreeNode(NodeType nodeType);
    RTreeNode(Triangle dataObject, BR dataBounds);
    ~RTreeNode();
};

class RTree
{
public:
    RTree();
    RTree(int MIN_CHILDREN, int MAX_CHILDREN);
    ~RTree();
    void Insert(Triangle &triangle);
    vector<Triangle> Search(pair<double, double> point) const; // Candidates for specific point

    // m is the minimum number of children allowed for each node
    // M is the maximum number of children allowed for each node.
    // Value of m is usually 40% of M and at most 50%.
    const int m_m = 2, m_M = 5;
    int m_depth = 1, m_dataNumber = 0;

private:
    vector<Triangle> Search(pair<double, double> &point, const RTreeNode &node) const;

    RTreeNode *ChooseLeaf(BR &triangleBounds, vector<RTreeNode *> &path) const;

    void SplitNode(RTreeNode *parentNode, RTreeNode *node);

    void PickSeeds(const RTreeNode *node, RTreeNode *&entry1, RTreeNode *&entry2) const;

    RTreeNode *PickNext(const vector<RTreeNode *> &entries, BR &bounds1, BR &bounds2, double &areaDiff) const;

    void AdjustTree(vector<RTreeNode *> &path, const BR &newBR);

    BR CalculateBounds(const BR &bounds1, const BR &bounds2) const;

    double CalculateBoundsArea(const BR &bounds) const;

    double AreaCost(const BR &bounds1, const BR &boundsOfInsertedNode) const;

    RTreeNode *root_; // Root node of the tree
};

ostream &operator<<(ostream &os, const Triangle &triangle);
ostream &operator<<(ostream &os, const BR &b);
ostream &operator<<(ostream &os, const RTree &tree);

#endif
