#ifndef NODE
#define NODE

#include <cmath>
#include "Point.h"

template <size_t N, class ElemType>
class Node {
public:
    // default constructor
    Node() {
        point = Point<N>();
        data = ElemType();
        left = nullptr;
        right = nullptr;
    }

    // fill constructor
    Node(const Point<N>& inpPoint, const ElemType& inpData) {
        point = inpPoint;
        data = inpData;
        left = nullptr;
        right = nullptr;
    }

    // copy constructor
    Node(const Node<N, ElemType>& constNode) {
        point = constNode.point;
        data = constNode.data;
        left = constNode.left;
        right = constNode.right;
    } 

    // pointers to the next two points in tree
    Node<N, ElemType>* left;
    Node<N, ElemType>* right;

    // associated point object
    Point<N> point;

    // associated value
    ElemType data;

};

template <size_t N, class ElemType>
bool operator==(Node<N, ElemType> node, Point<N> point) {
    for (int i = 0; i < point.size(); i++) {
        if (node.data[i] != point[i]) {
            return false;
        }
    }
    return true;
}

template <size_t N, class ElemType>
bool operator==(Point<N> point, Node<N, ElemType> node) {
    for (int i = 0; i < point.size(); i++) {
        if (node.data[i] != point[i]) {
            return false;
        }
    }
    return true;
}

template <size_t N, class ElemType>
std::ostream& operator<<(std::ostream& output, const Node<N, ElemType>& node) {
    output << "(";
    Point<N> point = node.point;
    for (double coord : point) {
        output << coord << ", ";
    }
    output << "\b\b)";
    return output;
}

#endif // NODE
