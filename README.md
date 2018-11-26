# KD Ttree
A kd-tree is a generalization of a binary search tree that stores points in k-dimensional space. That is, you could use a kd-tree to store a collection of points in the Cartesian plane, in three-dimensional space, etc. This program implements a k-dimensional binary search tree and also implements the k-nearest-neighbors search algorithm, which is an example of a **similar** search as opposed to the traditional **exact** search.

## Partition Dimension
The KD tree is implemented by having a different partition dimension, or coordinate index, at every level of the tree. For example, at the first level of the tree we'd only have one root node and we choose its children nodes by comparing the children nodes and the root nodes respective **first** coordinate values. At the second level, we choose children nodes for the third level by comparing the third level nodes with the second level nodes at the **second** coordinate value. This has the effect of splitting the n-dimensional space along a different dimension every time.

## KNN Search
The similarity search is implemented with a **bounded** priority queue (a queuu)with a comparator function that represents how close (Euclidean distance) to the target point the given point is. Points that are closer to the target point have a higher priority. The algorithm starts at the root node, enqueues it into the bounded priority queue, then decides whether or not to go into the left or right subtree.

It makes this decision by seeing which subtree the target point would be located in, according to the relevant partition dimension. It then repeats this process. It will not go into the other subtree unless either (1) the bounded priority queue is not yet full or (2) the distance from the target point to the furthest point in the initially chosen subtree is greater than the partition distance from the partition hyperplane to the point. In the second case, we couldn't be sure that these isn't a better point located in the other subtree. 

## Credits
This project was part of CS 106L (Standard C++ Programming) taught by Ali Malik in 2018 at Stanford University.
