#include "RTree.h"
#include <iostream>
#include <libspatialindex/libspatialindex.h>

using namespace std;

int main(int argc, char const *argv[])
{
    RTree tree(2);
    cout << tree.m_M << endl;
    cout << tree.children.trees.size() << endl;

    return 0;
}
