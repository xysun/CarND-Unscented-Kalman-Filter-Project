#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */


	if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	VectorXd rmse = VectorXd::Zero(estimations[0].size());
	
	for (int i = 0; i < estimations.size(); i++) {
		VectorXd diff = estimations[i] - ground_truth[i];
		diff = diff.array() * diff.array();
		rmse += diff;
	}

	// take mean
	rmse = rmse / (estimations.size());

	// sqrt
	rmse = rmse.array().sqrt();

	return rmse;
}
