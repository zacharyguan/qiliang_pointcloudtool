//
// Created by Zachary on 2023/10/12.
//

#ifndef PLYFASTLOADER_KDTREE_H
#define PLYFASTLOADER_KDTREE_H

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include "pct_entity.h"

namespace PCT {
    template<size_t DIM> using point_t = std::array< double, DIM >;
    using indexArr = std::vector< size_t >;
    template<size_t DIM> using pointIndex = typename std::pair< point_t<DIM>&, size_t >;
    template<size_t DIM> using pointIndexArr = typename std::vector< pointIndex<DIM> >;
    template<size_t DIM> using pointVec = std::vector< point_t<DIM> >;

    class KDNode {
    public:
        size_t leftIndex;
        size_t rightIndex;

        // initializer
        KDNode() {
            leftIndex = rightIndex = std::numeric_limits<size_t>::max();
        };
        ~KDNode() = default;
    };

    template<size_t DIM>
    class KDTree {
        std::vector<KDNode> tree;
        size_t root;
        const size_t leaf = std::numeric_limits<size_t>::max();
        IMultiDimensionDataSet& dataSet;
    public:
        explicit KDTree(IMultiDimensionDataSet &multiDimensionDataSet) : dataSet(multiDimensionDataSet) {
            if (DIM == 0) {
                throw std::runtime_error("DIM must large than 0.");
            }
            // iterators
            indexArr arr;
            arr.resize(dataSet.dataCount());
            std::iota(arr.begin(), arr.end(), 0);
            tree.resize(dataSet.dataCount());
            std::fill(tree.begin(), tree.end(), KDNode());

            auto begin = arr.begin();
            auto end = arr.end();

            size_t length = arr.size();
            size_t level = 0;  // starting

            root = KDTree::make_tree(begin, end, length, level);
        }
    private:
        size_t make_tree(const typename indexArr::iterator &begin,
                                 const typename indexArr::iterator &end,
                                 const size_t &length,
                                 const size_t &level) {
            if (begin == end) {
                return leaf;  // empty tree
            }

            size_t l_len = length / 2;
            size_t r_len = length - l_len - 1;

            auto middle = begin + l_len;

            if (length > 1) {
                std::nth_element(begin, middle, end, [this,&level](size_t& a, size_t& b) -> bool {
                    return (this->dataSet.dimensionValueAtIndex(a, level) < this->dataSet.dimensionValueAtIndex(b, level));
                });
            }

            size_t nextLevel = (level + 1) % DIM;

            size_t left;
            if (l_len > 0) {
                left = make_tree(begin, middle, l_len, nextLevel);
            } else {
                left = leaf;
            }
            size_t right;
            if (r_len > 0) {
                right = make_tree(middle + 1, end, r_len, nextLevel);
            } else {
                right = leaf;
            }

            return saveNode(*middle, left, right);
        }

        inline size_t saveNode(const size_t &index, const size_t &left, const size_t &right) {
            tree[index].leftIndex = left;
            tree[index].rightIndex = right;
            return index;
        }

        point_t<DIM> getPoint(const size_t index) const {
            point_t<DIM> p;
            for (int i=0; i<DIM; i++) {
                p[i] = dataSet.dimensionValueAtIndex(index, i);
            }
            return p;
        }

        // square euclidean distance
        inline double dist2(const size_t &aIndex, const point_t<DIM> &b) {
            double distc = 0;
            for (size_t i = 0; i < DIM; i++) {
                double di = dataSet.dimensionValueAtIndex(aIndex, i) - b.at(i);
                distc += di * di;
            }
            return distc;
        }

// euclidean distance
        size_t nearest_(
                const size_t &branch,
                const point_t<DIM> &pt,
                const size_t &level,
                const size_t &best,
                const double &best_dist) {
            double d, dx, dx2;

            if (branch == leaf) {
                return leaf;  // basically, null
            }

            size_t dim = DIM;

            d = dist2(branch, pt);
            dx = dataSet.dimensionValueAtIndex(branch, level) - pt.at(level);
            dx2 = dx * dx;

            size_t best_l = best;
            double best_dist_l = best_dist;

            if (d < best_dist) {
                best_dist_l = d;
                best_l = branch;
            }

            size_t next_lv = (level + 1) % dim;
            size_t section;
            size_t other;

            // select which branch makes sense to check
            if (dx > 0) {
                section = tree[branch].leftIndex;
                other = tree[branch].rightIndex;
            } else {
                section = tree[branch].rightIndex;
                other = tree[branch].leftIndex;
            }

            // keep nearest neighbor from further down the tree
            size_t further = nearest_(section, pt, next_lv, best_l, best_dist_l);
            if (further != leaf) {
                double dl = dist2(further, pt);
                if (dl < best_dist_l) {
                    best_dist_l = dl;
                    best_l = further;
                }
            }
            // only check the other branch if it makes sense to do so
            if (dx2 < best_dist_l) {
                further = nearest_(other, pt, next_lv, best_l, best_dist_l);
                if (further != leaf) {
                    double dl = dist2(further, pt);
                    if (dl < best_dist_l) {
                        best_dist_l = dl;
                        best_l = further;
                    }
                }
            }

            return best_l;
        }

        // default caller
        size_t nearest_(const point_t<DIM> &pt) {
            size_t level = 0;
            // KDNodePtr best = branch;
            double branch_dist = dist2(root, pt);
            return nearest_(root,          // beginning of tree
                            pt,            // point we are querying
                            level,         // start from level 0
                            root,          // best is the root
                            branch_dist);  // best_dist = branch_dist
        }
    public:
        point_t<DIM> nearest_point(const point_t<DIM> &pt) {
            return getPoint(nearest_(pt));
        }
        size_t nearest_index(const point_t<DIM> &pt) {
            return nearest_(pt);
        }
        pointIndex<DIM> nearest_pointIndex(const point_t<DIM> &pt) {
            size_t nearest = nearest_(pt);
            return pointIndex<DIM>(getPoint(nearest), nearest);
        }

    private:
        indexArr neighborhood_(
                const size_t &branch,
                const point_t<DIM> &pt,
                const double &rad,
                const size_t &level) {
            double d, dx, dx2;

            if (branch == leaf) {
                // branch has no point, means it is a leaf,
                // no points to add
                return {};
            }

            size_t dim = DIM;

            double r2 = rad * rad;

            d = dist2(branch, pt);
            dx = dataSet.dimensionValueAtIndex(branch, level) - pt.at(level);
            dx2 = dx * dx;

            indexArr nbh, nbh_s, nbh_o;
            if (d <= r2) {
                nbh.push_back(branch);
            }

            //
            size_t section;
            size_t other;
            if (dx > 0) {
                section = tree[branch].leftIndex;
                other = tree[branch].rightIndex;
            } else {
                section = tree[branch].rightIndex;
                other = tree[branch].leftIndex;
            }

            nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
            nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
            if (dx2 < r2) {
                nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
                nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
            }

            return nbh;
        }

        size_t neighborhood_dist_(
                const size_t &branch,
                const point_t<DIM> &pt,
                const double &rad,
                const size_t &level,
                std::vector<size_t>& indexes,
                std::vector<double>& distances) {
            double d2, dx, dx2;

            if (branch == leaf) {
                // branch has no point, means it is a leaf,
                // no points to add
                return 0;
            }

            size_t dim = DIM;

            double r2 = rad * rad;

            d2 = dist2(branch, pt);
            dx = dataSet.dimensionValueAtIndex(branch, level) - pt.at(level);
            dx2 = dx * dx;

            if (d2 <= r2) {
                indexes.emplace_back(branch);
                distances.emplace_back(std::sqrt(d2));
            }

            //
            size_t section;
            size_t other;
            if (dx > 0) {
                section = tree[branch].leftIndex;
                other = tree[branch].rightIndex;
            } else {
                section = tree[branch].rightIndex;
                other = tree[branch].leftIndex;
            }

            size_t nextLevel = (level + 1) % dim;
            neighborhood_dist_(section, pt, rad, nextLevel, indexes, distances);
            if (dx2 < r2) {
                neighborhood_dist_(other, pt, rad, nextLevel, indexes, distances);
            }

            return indexes.size();
        }

    public:
        pointIndexArr<DIM> neighborhood(
                const point_t<DIM> &pt,
                const double &rad) {
            size_t level = 0;
            indexArr nbh = neighborhood_(root, pt, rad, level);
            pointIndexArr<DIM> nbhp;
            nbhp.resize(nbh.size());
            std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                           [this](size_t index) { return pointIndexArr<DIM>(this->getPoint(index),index); });

            return nbhp;
        }

        size_t neighborhood(
                const point_t<DIM> &pt,
                const double &rad,
                std::vector<size_t>& indexes,
                std::vector<double>& distances) {
            size_t level = 0;
            return neighborhood_dist_(root, pt, rad, level, indexes, distances);
        }

        pointVec<DIM> neighborhood_points(
                const point_t<DIM> &pt,
                const double &rad) {
            size_t level = 0;
            indexArr nbh = neighborhood_(root, pt, rad, level);
            pointVec<DIM> nbhp;
            nbhp.resize(nbh.size());
            std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                           [this](size_t index) { return this->getPoint(index); });
            return nbhp;
        }

        indexArr neighborhood_indices(
                const point_t<DIM> &pt,
                const double &rad) {
            size_t level = 0;
            return neighborhood_(root, pt, rad, level);
        }
    };
}

#endif //PLYFASTLOADER_KDTREE_H
