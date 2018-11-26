/**
 * File: KDTree.h
 * Author: Kevin Tan
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "Node.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>

using namespace std;

template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();

    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;

private:

    // pointer to root node of KDTree
    Node<N, ElemType>* root;

    // helper functions
    size_t sizeHelper(Node<N, ElemType>* node) const;
    void insertHelper(const Node<N, ElemType>*& newPoint,
                      Node<N,ElemType>*& currNode, size_t depth);
    void destructorHelper(Node<N, ElemType>*& node);

    Node<N, ElemType>* findHelper(const Point<N>& pt, Node<N, ElemType>* curr, size_t depth) const;
    Node<N, ElemType>* find(const Point<N>& pt) const;

    void printKDTreeHelper(Node<N, ElemType>* node) const;
    void printKDTree() const;

    void kNNValueHelper(const Point<N>& key, Node<N, ElemType>* node, BoundedPQueue<Node<N, ElemType>>& bpq, int depth) const;

    void copyConstructorHelper(Node<N, ElemType>* node);

};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    root = nullptr;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::destructorHelper(Node<N, ElemType>*& node) {
    if (node != nullptr) {
        destructorHelper(node->left);
        destructorHelper(node->right);
        delete node;
    }
    // base case: root == nullptr; do nothing.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    destructorHelper(root);
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template<size_t N, class ElemType>
size_t KDTree<N, ElemType>::sizeHelper(Node<N, ElemType>* node) const {
    if (node == nullptr) {
        // base case
        return 0;
    } else {
        // recursive case
        return sizeHelper(node->left) + sizeHelper(node->right) + 1;
    }
}

template<size_t N, class ElemType>
size_t KDTree<N, ElemType>::size() const {
    return sizeHelper(root);
}

template<size_t N, class ElemType>
bool KDTree<N, ElemType>::empty() const {
    return root == nullptr;
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::insertHelper(const Node<N, ElemType>*& insertNode,
                                       Node<N, ElemType>*& currNode, size_t depth) {
    if (currNode == nullptr) {
        // if the KDTree is empty or you find an open space
        Node<N, ElemType>* newNode = new Node<N, ElemType>(*insertNode);
        currNode = newNode;
    } else if (currNode->point == insertNode->point) {
        // there is a repeat point: rewrite Elemtype
        currNode->data = insertNode->data;
    } else {
        // recursive case: continue searching tree for opening or repeat
        size_t partitionDim = depth % dimension();
        if (insertNode->point[partitionDim] <= currNode->point[partitionDim]) {
            insertHelper(insertNode, currNode->left, ++depth);
        } else {
            insertHelper(insertNode, currNode->right, ++depth);
        }
    }
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    const Node<N, ElemType>* newPoint = new Node<N, ElemType>(pt, value);
    insertHelper(newPoint, root, 0);
}

// Helper function for the find function.
template<size_t N, class ElemType>
Node<N, ElemType>* KDTree<N, ElemType>::findHelper(const Point<N>& pt, Node<N, ElemType>* curr, size_t depth) const {
    if (curr == nullptr || curr->point == pt) {
        // current node is empty or equal to point
        return curr;
    } else {
        // move on next level
        size_t partitionDimension = depth % pt.size();
        if (pt[partitionDimension] <= curr->point[partitionDimension]) {
            return findHelper(pt, curr->left, ++depth);
        } else {
            return findHelper(pt, curr->right, ++depth);
        }
    }
}

// Internal function that the following functions use internally:
//  1. contains
//  2. operator[]
//  3. at()
//  4. at() const
template<size_t N, class ElemType>
Node<N, ElemType>* KDTree<N, ElemType>::find(const Point<N>& pt) const {
    return findHelper(pt, root, 0);
}

template<size_t N, class ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    return (find(pt) != nullptr);
}

/*
 * Returns an Elemtype by reference located at a certain Point<N> object.
 * If the point doesn't exist, it is created.
 */
template<size_t N, class ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    Node<N, ElemType>* node = find(pt);
    if (node == nullptr) {
        // if node doesn't exist, put it in.
        insert(pt, ElemType());
        node = find(pt);
    }
    return node->data;
}

template<size_t N, class ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {

    // make current object into const object
    const KDTree<N, ElemType>& tree = *this;

    // call .at method on const object
    return const_cast<ElemType&>(tree.at(pt));

}

template<size_t N, class ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    Node<N, ElemType>* node = find(pt);
    if (node == nullptr) {
        // node not in tree
        throw out_of_range(" ");
    } else {
        // return value
        return const_cast<ElemType&>(node->data);
    }
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::printKDTreeHelper(Node<N, ElemType>* node) const {
    if (node == nullptr) {
        cout << "nullptr" << endl;
    } else {
        cout << "[" << node->point << " : " << node->data << "]" << endl;
        cout << "left: ";
        printKDTreeHelper(node->left);
        cout << "right: ";
        printKDTreeHelper(node->right);
    }
}

template<size_t N, class ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {
    BoundedPQueue<Node<N, ElemType>> bpq(k);

    // modify bpq to get closest neighbors
    int depth = 0;
    kNNValueHelper(key, root, bpq, depth);

    // create HashMap of frequencies
    map<ElemType, int> freq;

    // calculate frequency of labels
    while (!bpq.empty()) {
        Node<N, ElemType> node = bpq.dequeueMin();
        if (freq.count(node.data) == 0) {
            freq.emplace(node.data, 1);
        } else {
            freq[node.data]++;
        }
    }

    // return most frequently occuring label
    ElemType mostFrequentLabel = ElemType();
    int maxFreq = 0;
    for (auto elem : freq) {
        if (freq.at(elem.first) > maxFreq) {
            mostFrequentLabel = elem.first;
            maxFreq = freq.at(elem.first);
        }
    }

    return mostFrequentLabel;
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::kNNValueHelper(const Point<N>& key, Node<N, ElemType>* node, BoundedPQueue<Node<N, ElemType>>& bpq, int depth) const {
    if (node == nullptr) {
        return;
    }

    // add current node to BPQ
    bpq.enqueue(*node, Distance(key, node->point));

    // recursively search half of tree with test point
    int partitionDim = depth % key.size();
    if (key[partitionDim] < node->point[partitionDim]) {
        kNNValueHelper(key, node->left, bpq, depth++);
    } else {
        kNNValueHelper(key, node->right, bpq, depth++);
    }

    // if candidate hypersphere crosses splitting plane
    // or bpq is not full, explore the other branch
    if (bpq.size() < bpq.maxSize() || abs(node->point[partitionDim] - key[partitionDim]) < bpq.worst()) {
        if (key[partitionDim] < node->point[partitionDim]) {
            kNNValueHelper(key, node->right, bpq, depth++);
        } else {
            kNNValueHelper(key, node->left, bpq, depth++);
        }
    }
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::copyConstructorHelper(Node<N, ElemType>* node) {
    if (node == nullptr) {
        return;
    }
    node = new Node<N, ElemType>(*node);
    copyConstructorHelper(node->left);
    copyConstructorHelper(node->right);
}

// performs deep copy into a new KDTree object
// note that you are allowed to access private variables - sibling access
template<size_t N, class ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    root = new Node<N, ElemType>(*(rhs.root));
    copyConstructorHelper(root->left);
    copyConstructorHelper(root->right);
}

// performs a deep copy into an old KDTree object
// (1) clear up previous data
// (2) make sure to handle self-reference
// (3) return a reference to the object
template<size_t N, class ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
    if (this != &rhs) {
        // clear up old data
        destructorHelper(root);
        root = new Node<N, ElemType>(*(rhs.root));
        copyConstructorHelper(root->left);
        copyConstructorHelper(root->right);
    }
    // else, self-reference: do nothing
    return *this;
}

template<size_t N, class ElemType>
void KDTree<N, ElemType>::printKDTree() const {
    printKDTreeHelper(root);
}

#endif // KDTREE_INCLUDED
