#include <cmath>
#include "RTree.h"

using namespace std;

Triangle::Triangle() {}

Triangle::Triangle(const double x0, const double y0, double x1, double y1, double x2, double y2) : vertex0(make_pair(x0, y0)), vertex1(make_pair(x1, y1)), vertex2(make_pair(x2, y2)) {}

bool Triangle::PointInTriangle(pair<double, double> &point) const
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

ostream &operator<<(ostream &os, const RTree &tree)
{
    return os << "L'arbre de paramètres (m, M) = (" << tree.m_m << ", " << tree.m_M << ") contient " << tree.m_dataNumber << " objets et a une hauteur de " << tree.m_depth << endl;
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
    for (auto &child : this->children)
        delete child;
}

BR MBR_of_triangle(Triangle &triangle) // Minimum Bounding Rectangle
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

vector<Triangle> RTree::Search(pair<double, double> &point, const RTreeNode &node) const
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

void RTree::Insert(Triangle &triangle)
{
    BR triangleBounds = MBR_of_triangle(triangle);
    // Find the leaf node where the point should be inserted
    vector<RTreeNode *> path;
    RTreeNode *leaf = ChooseLeaf(triangleBounds, path);

    // Add the triangle to the leaf node's children
    RTreeNode *newNode = new RTreeNode(triangle, triangleBounds);
    leaf->children.push_back(newNode);
    m_dataNumber++;
    // Update the bounding box for the leaf node to include the new point
    AdjustTree(path, triangleBounds);
}

// Returns a vector starting from the root to the leaf
RTreeNode *RTree::ChooseLeaf(BR &triangleBounds, vector<RTreeNode *> &path) const
{
    RTreeNode *node = root_; // Start at the root of the tree
    path.push_back(node); // We keep track of each node visited starting with the root node

    while (node->nodeType == INTERNAL)
    {
        // Choose the child node whose rectangle needs the least enlargement to include triangle
        double min_areaCost = INFINITY;
        // Child candidate
        RTreeNode *min_child = node->children[0];
        for (const auto &child : node->children)
        {
            // Enlargement needed to include triangle
            double areaCost = AreaCost(child->bounds, triangleBounds);
            // Resovle ties by choosing the child with the rectangle of smallest area
            if (areaCost < min_areaCost || (areaCost == min_areaCost && CalculateBoundsArea(child->bounds) < CalculateBoundsArea(min_child->bounds)))
            {
                min_areaCost = areaCost;
                min_child = child;
            }
        }
        // Chosen node
        node = min_child;
        path.push_back(node);
    }
    return node;
}

void RTree::SplitNode(RTreeNode *parentNode, RTreeNode *node)
{
    // Choose the two entries to be the first elements of the groups
    RTreeNode *seed1, *seed2;
    PickSeeds(node, seed1, seed2);
    // Create two new leaf nodes in which the entries will be distributed
    RTreeNode *node1 = new RTreeNode(node->nodeType);
    RTreeNode *node2 = new RTreeNode(node->nodeType);
    // Add both seeds to one of the new leafs, update leafs' bounds and remove the seeds from the remaining entries
    node1->children = {seed1};
    node1->bounds = seed1->bounds;
    node->children.erase(remove(node->children.begin(), node->children.end(), seed1), node->children.end());

    node2->children = {seed2};
    node2->bounds = seed2->bounds;
    node->children.erase(remove(node->children.begin(), node->children.end(), seed2), node->children.end());

    while (!node->children.empty())
    {   
        // If one group has so few entries that all the rest must be assigned to it in order for it to have the minimum number m_m, assign them and stop
        if (node->children.size() + node1->children.size() <= m_m)
        {
            for (const auto &child : node->children)
            {
                node1->children.push_back(child);
                node1->bounds = CalculateBounds(node1->bounds, child->bounds);
            }
            node->children = {};
            break;
        }

        if (node->children.size() + node2->children.size() <= m_m)
        {
            for (const auto &child : node->children)
            {
                node2->children.push_back(child);
                node2->bounds = CalculateBounds(node2->bounds, child->bounds);
            }
            node->children = {};
            break;
        }

        // Choose next entry to assign
        double areaDiff;
        RTreeNode *entry = PickNext(node->children, node1->bounds, node2->bounds, areaDiff);

        double area1 = CalculateBoundsArea(node1->bounds),
               area2 = CalculateBoundsArea(node2->bounds);
        if (areaDiff > 0 || (areaDiff == 0 && area1 < area2) || (area1 == area2 && node1->children.size() <= node2->children.size()))
        {
            node1->children.push_back(entry);
            node1->bounds = CalculateBounds(node1->bounds, entry->bounds);
        }
        else
        {
            node2->children.push_back(entry);
            node2->bounds = CalculateBounds(node2->bounds, entry->bounds);
        }
        node->children.erase(remove(node->children.begin(), node->children.end(), entry), node->children.end());
    }

    // Make the original node an internal node and add the two new nodes as its children
    // node->nodeType = INTERNAL;
    parentNode->children.erase(remove(parentNode->children.begin(), parentNode->children.end(), node), parentNode->children.end());
    delete node;
    parentNode->children.push_back(node1);
    parentNode->children.push_back(node2);
}

void RTree::PickSeeds(const RTreeNode *node, RTreeNode *&entry1, RTreeNode *&entry2) const
{
    // Along each dimension, find the entry whose rectangle has the highest low side, and the one with the lowest high side.
    RTreeNode *resultLeftNode;
    RTreeNode *resultRightNode;
    double normalizedSep = -INFINITY;
    for (int i = 0; i < node->bounds.size(); i++)
    {
        entry1 = node->children[0];
        for (const auto &child : node->children)
        {
            if (child->bounds[i].second < entry1->bounds[i].second)
            {
                entry1 = child;
            }
        }

        entry2 = entry1 != node->children[0] ? node->children[0] : node->children[1];
        for (const auto &child : node->children)
        {
            if (child == entry1)
                continue;

            if (child->bounds[i].first > entry2->bounds[i].first)
            {
                entry2 = child;
            }
        }

        // Normalize the separations by dividing by the width of the entire set along the corresponding dimension
        double L = max(entry1->bounds[i].first - entry2->bounds[i].second, entry2->bounds[i].first - entry1->bounds[i].second);
        double W = node->bounds[i].second - node->bounds[i].first;
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

RTreeNode *RTree::PickNext(const vector<RTreeNode *> &entries, BR &bounds1, BR &bounds2, double &areaDiff) const
{
    double max_areaIncreaseDiff = -INFINITY;
    RTreeNode *max_entry;
    for (const auto &entry : entries)
    {
        double d1 = AreaCost(bounds1, entry->bounds);
        double d2 = AreaCost(bounds2, entry->bounds);
        if (abs(d2 - d1) > max_areaIncreaseDiff)
        {
            max_areaIncreaseDiff = abs(d2 - d1);
            areaDiff = d2 - d1;
            max_entry = entry;
        }
    }

    return max_entry;
}

void RTree::AdjustTree(vector<RTreeNode *> &path, const BR &newBR)
{
    // From leaf to root
    while (!path.empty())
    {
        RTreeNode *node = path.back();
        path.pop_back();

        // If the node has too many children, split it into two nodes
        if (node->children.size() > m_M)
        {
            if (!path.empty())
            {
                SplitNode(path.back(), node);
            }
            else
            {
                RTreeNode *newRoot = new RTreeNode(INTERNAL);
                SplitNode(newRoot, node);
                root_ = newRoot;
                path.push_back(root_);
                m_depth++;
            }
            // node's bounds were already updated by SplitNode
            continue;
        }

        node->bounds = CalculateBounds(node->bounds, newBR);
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
    // Initial area
    double area = CalculateBoundsArea(bounds1);

    // Area when adding new node
    BR newBounds = CalculateBounds(bounds1, boundsOfInsertedNode);
    // Difference between new area and initial area
    return CalculateBoundsArea(newBounds) - area;
}