#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
namespace statistic{
    std::vector<double> logistic_regression
        (const std::vector< std::vector<double> >& data, const std::vector<bool> is_treated);
    double predict(const std::vector<double> & weights, const std::vector<double> & property);

    double cross_validation_5_fold(std::vector< std::vector<double> >& data, const std::vector<bool> gound_truth);

    std::vector<double> ada_boosting(std::vector< std::vector<double> >& data, std::vector< std::vector<double> >& learner, const std::vector<bool> gound_truth);
}