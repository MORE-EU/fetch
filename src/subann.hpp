#include <unordered_map>
#include <string>
#include "utils.hpp"
#include "timeseries.hpp"
#include <omp.h>
#include <pybind11/pybind11.h>
#include <map>
namespace py = pybind11;
class SubAnnDS {
private:
    TimeSeries signal;
    double w1, s1, w2, s2, delta, window, leq;
    //double w1, s1, delta;
    std::unordered_map<std::vector<int>, std::vector<std::pair<int,int>>, VectorHasher> dict;
    //std::map<std::pair<int, int>, bool> candidates;

public:
    SubAnnDS(TimeSeries &signal, const double w1, const double delta);
    SubAnnDS(TimeSeries& data, const double w1, const double delta, const double window, const bool leq);
    SubAnnDS(TimeSeries& data, const double w1, const double w2, const double delta, const double window, const bool leq);

    const std::unordered_map<std::vector<int>, std::vector<std::pair<int,int>>, VectorHasher> get_dict() const {
        return dict;
    }
    const std::vector<std::pair<int,int>> getAnswersIds1D(TimeSeries query);
    std::vector<std::pair<int, int>> range_query(TimeSeries query, double R);
    std::vector<std::pair<int, int>> exact_range_query(TimeSeries query, double R);

    std::pair<int, int>  ann_query(TimeSeries query);
    std::pair<int, int>  exact_nn_query(TimeSeries query);

    friend ostream& operator<<(ostream& os, const SubAnnDS& obj);
    void stats();
};