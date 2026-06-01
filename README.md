# Maximize total distance given order and radius 
Given integers $n$ and $r$, this repository contains code to calculate the maximum total distance among
(1) All $n$-vertex graphs of radius $r$. The extremal graph is a tree and the corresponding code can be found in "uniqueCenter". This code also computes the extremal graph among all $n$-vertex trees of diameter $2r$.
(2) All $n$-vertex trees of diameter $2r-1$. The extremal graph is a tree whose center contains exactly two vertices and the corresponding code can be found in "twoCenters".
(3) All $n$-vertex clique-rooted trees of radius $r$. The corresponding code can be found in "cliqueCenter".


The program can be compiled by navigating to the directory and then:
```bash
g++ -g -std=c++11 -O3 computeMaxWienerIndexGivenOrderAndRadius.cpp -o computeMaxWienerIndexGivenOrderAndRadiusExecutable
```

An example of how the program can be used:
```bash
cd uniqueCenter/
./computeMaxWienerIndexGivenOrderAndRadiusExecutable 100 6
```
Computes the maximum total distance among all trees on 100 vertices that have radius 6. The output of the program is as follows:

```bash
Start output for case where the tree has a single centre:
Maximum Wiener index of a graph with 100 vertices and radius 6 is: 43980
Take dW(100,19,6) 1 times
Leave counts: 6 7 
New vertex with labels 0 1 at distance 2 from the leaves

Take dW(100,20,6) 4 times
Leave counts: 7 7 
New vertex with labels 0 1 at distance 2 from the leaves

Finished output for case where the tree has a single centre.
```

This means that the maximum total distance is equal to 43980. 

The optimal trees is also displayed in a succinct manner. In this case, the root of the optimal tree has degree 5. One branch has 19 vertices, whereas the other four branches each have 20 vertices. The branches themselves are also displayed in a succinct manner by displaying the number of leaves and displaying where the vertices of degree 3 are located. Vertices of degree 2 are not displayed to ensure this output remains short and interpretable.
