// Independent exact skeleton DP. C++20 and GMP; all child pairs are tested.
// Unlike the supplied description, this does NOT use an interval-overlap sweep.
// Build: g++ -O3 -std=c++20 skeleton_envelope.cpp -lgmpxx -lgmp -o skeleton_envelope
// Run: ./skeleton_envelope 101
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using Z = mpz_class;
using Q = mpq_class;
constexpr size_t NONE = std::numeric_limits<size_t>::max();
struct Node {
    Z b;
    Q R;
    size_t left = NONE, right = NONE, shape = 0;
    bool symmetric = true;
};
struct Candidate { Z b; Q R; size_t left, right; Q start; };
std::vector<Node> nodes;
std::map<std::pair<size_t,size_t>,size_t> shape_ids;

// Prove a common finite support window. For q=r-1, put
// C <= (sum_{k=1}^{q-1} sqrt(k))^2, using dyadic lower bounds.
// Cauchy-Schwarz gives R >= 1+C/(b-q). A feasible single-fork
// skeleton gives an upper bound U on mu. Hence any minimizer satisfies
// b^2+(C-q-U)b+Uq <= 0. Find its largest possible integer b, B.
// Its supporting slope y=b/R is at most B/(1+C/(B-q)).
Q support_window(int max_r) {
    Q window = 1;
    const Z scale = Z(1) << 32;
    for (long q=2; q<max_r; ++q) {
        Z sum = 0;
        for (long k=1; k<q; ++k) {
            Z arg = Z(k)*scale*scale, root;
            mpz_sqrt(root.get_mpz_t(), arg.get_mpz_t());
            sum += root;
        }
        Q C(sum*sum, scale*scale); C.canonicalize();
        Q U; Z witness_b; bool first=true;
        for (long k=1; k<=q; ++k) {
            // k common edges, then two arms of length q-k.
            Z b = Z(k)*(q+1)-Z(k)*(k+1)/2+Z(q-k)*(q-k+1);
            Q R(q+k,2); R.canonicalize();
            Q product = Q(b)*R;
            if (first || product<U) { U=product; witness_b=b; first=false; }
        }
        auto polynomial = [&](const Z& b)->Q {
            return Q(b*b)+(C-q-U)*Q(b)+U*q;
        };
        if (polynomial(witness_b)>0) throw std::runtime_error("bad bound");
        Z lo=witness_b, hi=U.get_num()/U.get_den()+2;
        while (hi-lo>1) {
            Z middle=(lo+hi)/2;
            if (polynomial(middle)<=0) lo=middle; else hi=middle;
        }
        Q R_lower=1+C/Q(lo-q);
        Q bound=Q(lo)/R_lower;
        if (bound>window) window=bound;
    }
    // An integer upper bound keeps subsequent comparisons simpler.
    return Q(window.get_num()/window.get_den()+1);
}

std::vector<size_t> next_level(const std::vector<size_t>& old,
                               long q, const Q& window) {
    std::map<Z,Candidate> by_b;
    auto insert = [&](const Z& b,const Q& R,size_t left,size_t right) {
        auto it=by_b.find(b);
        if (it==by_b.end()) by_b.emplace(b,Candidate{b,R,left,right,0});
        else if (R<it->second.R) it->second=Candidate{b,R,left,right,0};
    };
    for (size_t i:old) insert(nodes[i].b+q,nodes[i].R+1,i,NONE);
    for (size_t a=0;a<old.size();++a) {
        const auto &x=nodes[old[a]];
        for (size_t b=a;b<old.size();++b) {
            const auto &y=nodes[old[b]];
            insert(x.b+y.b+q,1+x.R*y.R/(x.R+y.R),old[a],old[b]);
        }
    }
    // Lower envelope of lines b+y*R, restricted to 0 <= y <= window.
    // States are considered by increasing b. Nondominated R decrease.
    std::vector<Candidate> hull;
    for (const auto &[b,candidate]:by_b) {
        Candidate c=candidate;
        if (!hull.empty() && c.R>=hull.back().R) continue;
        while (!hull.empty()) {
            c.start=Q(c.b-hull.back().b)/(hull.back().R-c.R);
            if (c.start>hull.back().start) break;
            hull.pop_back();
        }
        if (hull.empty()) c.start=0;
        if (c.start<=window) hull.push_back(c);
    }
    std::vector<size_t> result;
    for (const auto &c:hull) {
        size_t a=nodes[c.left].shape;
        size_t b=c.right==NONE ? NONE : nodes[c.right].shape;
        bool symmetric=nodes[c.left].symmetric;
        if (c.right!=NONE) {
            symmetric=symmetric && nodes[c.right].symmetric && a==b;
            if (b<a) std::swap(a,b);
        }
        auto [it,unused]=shape_ids.emplace(std::make_pair(a,b),shape_ids.size()+1);
        result.push_back(nodes.size());
        nodes.push_back(Node{c.b,c.R,c.left,c.right,it->second,symmetric});
    }
    return result;
}

// Independently traverse the returned tree and recompute b and R.
std::pair<Z,Q> check_shape(size_t id,long q) {
    const auto &v=nodes[id];
    if (v.left==NONE) {
        if (q!=1) throw std::runtime_error("terminal at wrong depth");
        return {Z(1),Q(1)};
    }
    auto a=check_shape(v.left,q-1);
    if (v.right==NONE) return {a.first+q,a.second+1};
    auto b=check_shape(v.right,q-1);
    return {a.first+b.first+q,1+a.second*b.second/(a.second+b.second)};
}

int main(int argc,char **argv) {
    try {
        int max_r=argc>1 ? std::stoi(argv[1]) : 101;
        if (max_r<2 || max_r>200) throw std::runtime_error("radius must be 2..200");
        Q window=support_window(max_r);
        std::cerr << "Certified support window: 0 <= y <= " << window << '\n';
        nodes.push_back(Node{Z(1),Q(1)});
        std::vector<size_t> states={0};
        std::cout << "r\tenvelope_states\tb\tR\tmu\tsymmetric\n";
        for (int r=2;r<=max_r;++r) {
            if (r>2) states=next_level(states,r-1,window);
            Q best=Q(nodes[states[0]].b)*nodes[states[0]].R;
            for (size_t id:states) {
                Q value=Q(nodes[id].b)*nodes[id].R;
                if (value<best) best=value;
            }
            for (size_t id:states) {
                const auto &v=nodes[id];
                if (Q(v.b)*v.R!=best) continue;
                auto checked=check_shape(id,r-1);
                if (checked.first!=v.b || checked.second!=v.R)
                    throw std::runtime_error("backtracking check failed");
                std::cout << r << '\t' << states.size() << '\t' << v.b << '\t'
                          << v.R << '\t' << best << '\t'
                          << (v.symmetric?"yes":"no") << std::endl;
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n'; return 1;
    }
}
