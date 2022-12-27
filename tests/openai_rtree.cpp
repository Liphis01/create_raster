class RTree
{
public:
    RTree() : root_(nullptr) {}

    // insert a new data point into the R-tree
    void insert(const Point &p, int id)
    {
        if (!root_)
        {
            // If the tree is empty, create a new leaf node as the root
            root_ = new Node(true);
        }

        // Choose the subtree that will contain the data point
        Node *n = root_->choose_subtree(p, 0);
        insert(p, id, n, 1);
    }

    // search for data points within a given rectangle
    std::vector<int> search(const Rect &r) const
    {
        std::vector<int> results;
        if (root_)
        {
            search(r, root_, results);
        }
        return results;
    }

private:
    // Recursive helper function for inserting a data point
    void insert(const Point &p, int id, Node *n, int level)
    {
        if (n->is_leaf())
        {
            // Add the data point to this leaf node
            n->add_child(p, id);
            return;
        }

        // Choose the subtree that will contain the data point
        Node *next = n->choose_subtree(p, level);
        insert(p, id, next, level + 1);
    }

    // Recursive helper function for searching for data points within a rectangle
    void search(const Rect &r, const Node *n, std::vector<int> &results) const
    {
        if (n->is_leaf())
        {
            // If this is a leaf node, check if its data points are within the rectangle
            for (const auto &child : n->children())
            {
                if (r.contains(child.point()))
                {
                    results.push_back(child.id());
                }
            }
        }
        else
        {
            // If this is an internal node, search its child nodes
            for (const auto &child : n->children())
            {
                if (r.intersects(child.bounds()))
                {
                    search(r, child, results);
                }
            }
        }
    }

    Node *root_; // Pointer to the root node of the R-tree
};