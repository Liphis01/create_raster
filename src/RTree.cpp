#include <cmath>
#include "RTree.h"

using namespace std;

Triangle::Triangle() {}

Triangle::Triangle(const double x0, const double y0, double x1, double y1, double x2, double y2) : vertex0(make_pair(x0, y0)), vertex1(make_pair(x1, y1)), vertex2(make_pair(x2, y2)) {}

bool Triangle::PointInTriangle(pair<double, double> point) const
{
    // Checks if a point is within a triangle's bounds
    double tx0 = this->vertex0.first, ty0 = this->vertex0.second;
    double tx1 = this->vertex1.first, ty1 = this->vertex1.second;
    double tx2 = this->vertex2.first, ty2 = this->vertex2.second;
    double px = point.first, py = point.second;

    double angle0, angle1, angle2;
    angle0 = (tx0 - px) * (ty1 - py) - (ty0 - py) * (tx1 - px);
    angle1 = (tx1 - px) * (ty2 - py) - (ty1 - py) * (tx2 - px);
    angle2 = (tx2 - px) * (ty0 - py) - (ty2 - py) * (tx0 - px);
    return (angle0 * angle1) >= 0 && (angle1 * angle2) >= 0;
}

double Triangle::Area() const
{
    return .5 * abs(vertex0.first * (vertex1.second - vertex2.second) + vertex1.first * (vertex2.second - vertex0.second) + vertex2.first * (vertex0.second - vertex1.second));
}

ostream &operator<<(ostream &os, const Triangle &triangle)
{
    return os << "Vertex0: (" << triangle.vertex0.first
              << ", " << triangle.vertex0.second << ")" << endl
              << "Vertex1: (" << triangle.vertex1.first
              << ", " << triangle.vertex1.second << ")" << endl
              << "Vertex2: (" << triangle.vertex2.first
              << ", " << triangle.vertex2.second << ")" << endl;
}

ostream &operator<<(ostream &os, const BR &b)
{
    for (int i = 0; i < b.size(); i++)
    {
        os << "Interval in dimension " << i << ": [" << b[i].first << ", " << b[i].second << "]" << endl;
    }
    return os;
}

RTreeNode::RTreeNode(NodeType nodeType) : nodeType(nodeType) {}

RTreeNode::RTreeNode(Triangle dataObject, BR dataBounds) : nodeType(DATA), dataObject(dataObject), bounds(dataBounds) {}

bool PointInBounds(pair<double, double> point, BR bounds)
{
    return bounds[0].first <= point.first && point.first <= bounds[0].second && bounds[1].first <= point.second && point.second <= bounds[1].second;
}

RTreeNode::~RTreeNode()
{
    // Recursively frees memory in the tree
    if (!this->children.empty())
        for (auto &child : this->children)
            delete child;
}

BR MBR_of_triangle(Triangle triangle) // Minimum Bounding Rectangle
{
    BR result;

    double min_x = min(min(triangle.vertex0.first, triangle.vertex1.first), triangle.vertex2.first);
    double max_x = max(max(triangle.vertex0.first, triangle.vertex1.first), triangle.vertex2.first);

    double min_y = min(min(triangle.vertex0.second, triangle.vertex1.second), triangle.vertex2.second);
    double max_y = max(max(triangle.vertex0.second, triangle.vertex1.second), triangle.vertex2.second);

    result.push_back(make_pair(min_x, max_x));
    result.push_back(make_pair(min_y, max_y));

    return result;
}

RTree::RTree()
{
    // Root node of the tree
    root_ = new RTreeNode(LEAF);
}

RTree::RTree(int MIN_CHILDREN, int MAX_CHILDREN) : m_m(MIN_CHILDREN), m_M(MAX_CHILDREN)
{
    // Root node of the tree
    root_ = new RTreeNode(LEAF);
}

RTree::~RTree()
{
    // Calls RTreeNode's destructor to recursively free memory allocated to the tree
    delete root_;
}

vector<Triangle> RTree::Search(pair<double, double> point) const
{
    return this->Search(point, *root_);
}

vector<Triangle> RTree::Search(pair<double, double> point, const RTreeNode &node) const
{
    if (node.nodeType == DATA && node.dataObject.PointInTriangle(point))
    {
        return {node.dataObject};
    }

    vector<Triangle> result = {};

    for (const auto &child : node.children)
    {
        if (PointInBounds(point, child->bounds))
        {
            vector<Triangle> v = this->Search(point, *child);
            // If a triangle is already found we stop searching (if we're interested in searching all triangles then remove this if statement and uncomment the line right after)
            if (!v.empty())
                return v;
            // result.insert(result.end(), v.begin(), v.end()); // Concatenation of v into result
        }
    }

    return result;
}

void RTree::Insert(Triangle triangle)
{
    BR triangleBounds = MBR_of_triangle(triangle);
    // Find the leaf node where the point should be inserted
    vector<RTreeNode *> path;
    RTreeNode *leaf = FindLeaf(triangleBounds, path);

    // Add the triangle to the leaf node's children
    RTreeNode *newNode = new RTreeNode(triangle, triangleBounds);
    leaf->children.push_back(newNode);
    // If the leaf node has too many children, split it into two nodes
    if (leaf->children.size() > m_M)
    {
        SplitNode(leaf);
    }
    // Update the bounding box for the leaf node to include the new point
    AdjustTree(path);
}

RTreeNode *RTree::FindLeaf(BR triangleBounds, vector<RTreeNode *> &path) const
{
    RTreeNode *node = root_; // Start at the root of the tree
    path.push_back(node); // We keep track of each node visited starting with the root node

    while (node->nodeType == INTERNAL)
    {
        // Choose the child node with the smallest cost in bounds
        double min_areaCost = INFINITY;
        RTreeNode *min_child = nullptr;
        for (const auto &child : node->children)
        {
            double areaCost = AreaCost(child->bounds, triangleBounds);
            // Resovle ties by choosing the child with the rectangle of smallest area
            if (areaCost < min_areaCost || (areaCost == min_areaCost && CalculateBoundsArea(child->bounds) < CalculateBoundsArea(min_child->bounds)))
            {
                min_areaCost = areaCost;
                min_child = child;
            }
        }
        node = min_child;
        path.push_back(node);
    }
    return node;
}

void RTree::SplitNode(RTreeNode *node)
{
    // Choose the two seeds that will be used to divide the entries into two groups
    RTreeNode *seed1, *seed2;
    ChooseSplitEntries(node, seed1, seed2);
    // Create two new leaf nodes
    RTreeNode *node1 = new RTreeNode(LEAF);
    RTreeNode *node2 = new RTreeNode(LEAF);
    // Divide the children of the original node between the two new nodes
    node1->children = {seed1};
    node1->bounds = seed1->bounds;
    node->children.erase(remove(node->children.begin(), node->children.end(), seed1), node->children.end());

    node2->children = {seed2};
    node2->bounds = seed2->bounds;
    node->children.erase(remove(node->children.begin(), node->children.end(), seed2), node->children.end());

    // Sort based on difference of area cost between adding to group 1 and group 2
    std::sort(node->children.begin(), node->children.end(), [&](RTreeNode *x, RTreeNode *y)
              { return AreaCost(node1->bounds, x->bounds) - AreaCost(node2->bounds, x->bounds) < AreaCost(node1->bounds, y->bounds) - AreaCost(node2->bounds, y->bounds); });
    while (node1->children.size() < m_m)
    {
        node1->children.push_back(node->children[0]);
        node1->bounds = CalculateBounds(node1->bounds, node->children[0]->bounds);
        node->children.erase(node->children.begin());
    }

    // Sort based on difference of area cost between adding to group 1 and group 2
    std::sort(node->children.begin(), node->children.end(), [&](RTreeNode *x, RTreeNode *y)
              { return AreaCost(node2->bounds, x->bounds) - AreaCost(node1->bounds, x->bounds) < AreaCost(node2->bounds, y->bounds) - AreaCost(node1->bounds, y->bounds); });
    while (node2->children.size() < m_m)
    {
        node2->children.push_back(node->children[0]);
        node2->bounds = CalculateBounds(node2->bounds, node->children[0]->bounds);
        node->children.erase(node->children.begin());
    }

    for (const auto &child : node->children)
    {
        // Calculate the difference between the new area (when adding child to each group) to the initial area
        double d1 = AreaCost(node1->bounds, child->bounds);
        double d2 = AreaCost(node2->bounds, child->bounds);

        // Add child to the node that costs the least area
        // Resolve ties by adding the seed to the group with fewer entries, then to either
        if (d1 < d2 || (d1 == d2 && (node1->children.size() <= node2->children.size())))
        {
            node1->children.push_back(child);
            node1->bounds = CalculateBounds(node1->bounds, child->bounds);
        }
        else
        {
            node2->children.push_back(child);
            node2->bounds = CalculateBounds(node2->bounds, child->bounds);
        }
    }
    // Make the original node an internal node and add the two new nodes as its children
    node->nodeType = INTERNAL;
    node->children = {node1, node2};
    node->bounds = CalculateBounds(node1->bounds, node2->bounds);
}

void RTree::ChooseSplitEntries(const RTreeNode *leaf, RTreeNode *&entry1, RTreeNode *&entry2) const
{
    // Along each dimension, find the entry whose rectangle has the highest low side, and the one with the lowest high side.
    RTreeNode *resultLeftNode = leaf->children[0];
    RTreeNode *resultRightNode = leaf->children[0];
    double normalizedSep = -INFINITY;
    for (int i = 0; i < leaf->bounds.size(); i++)
    {
        double min_highSide = INFINITY;
        double max_lowSide = -INFINITY;
        for (const auto &child : leaf->children)
        {
            if (child->bounds[i].second < min_highSide)
            {
                min_highSide = child->bounds[i].second;
                entry1 = child;
            }
            if (child->bounds[i].first > max_lowSide || (entry1 == entry2 && child->bounds[i].first == max_lowSide))
            {
                max_lowSide = child->bounds[i].first;
                entry2 = child;
            }
        }
        // Normalize the separations by dividing by the width of the entire set along the corresponding dimension
        double L = max_lowSide - min_highSide;
        double W = leaf->bounds[i].second - leaf->bounds[i].first;
        if (L / W > normalizedSep)
        {
            normalizedSep = L / W;
            resultLeftNode = entry1;
            resultRightNode = entry2;
        }
    }
    // Choose the pair with the greatest normalized separation along any dimension in the node's children
    entry1 = resultLeftNode;
    entry2 = resultRightNode;
}

void RTree::AdjustTree(vector<RTreeNode *> &path)
{
    while (!path.empty())
    {
        RTreeNode *node = path.back();
        path.pop_back();

        node->bounds = {};
        for (const auto &child : node->children)
        {
            node->bounds = CalculateBounds(node->bounds, child->bounds);
        }
    }
}

// Merge two bounds together
BR RTree::CalculateBounds(const BR &bounds1, const BR &bounds2) const
{
    if (bounds1.size() == 0)
        return bounds2;
    if (bounds2.size() == 0)
        return bounds1;

    double min_x = min(bounds1[0].first, bounds2[0].first);
    double max_x = max(bounds1[0].second, bounds2[0].second);
    double min_y = min(bounds1[1].first, bounds2[1].first);
    double max_y = max(bounds1[1].second, bounds2[1].second);
    return {make_pair(min_x, max_x), make_pair(min_y, max_y)};
}

double RTree::CalculateBoundsArea(const BR &bounds) const
{
    double width = bounds[0].second - bounds[0].first;
    double height = bounds[1].second - bounds[1].first;
    return width * height;
}

double RTree::AreaCost(const BR &bounds1, const BR &boundsOfInsertedNode) const
{
    double area = CalculateBoundsArea(bounds1);

    // Calculate the difference between the new area (when adding child to each group) to the initial area
    BR newBounds = CalculateBounds(bounds1, boundsOfInsertedNode);
    return CalculateBoundsArea(newBounds) - area;
}