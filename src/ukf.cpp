#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  
  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  time_us_ = 0;

  weights_ = VectorXd(2 * n_aug_ + 1);
  //set weights
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);


  NIS_laser_ = 0;
  NIS_radar_ = 0;

  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0,0,
	  0, 1, 0, 0,0;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
	  0, std_laspy_ * std_laspy_;

  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

	// initialize x_ and P_ with first measurement
	if (!is_initialized_) {

		// initiate x_ with first measurement;
		time_us_ = meas_package.timestamp_;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			double rho = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			double drho = meas_package.raw_measurements_[2];

			double px = rho * cos(phi);
			double py = rho * sin(phi);

			// initiate turn and turn rate to 0

			x_ << px, py, drho, phi, 0;

		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			laser only detect px,py
			*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}

		// initialize covariance
		P_ << 1, 0, 0, 0,0,
			0, 1, 0, 0,0,
			0, 0, 1000, 0,0,
			0, 0, 0, 1, 0,
			0, 0, 0, 0, 1;


		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	// get delta_t
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	// predict
	Prediction(delta_t);

	// update
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else {
		// laser
		UpdateLidar(meas_package);
	}

	// print the output
	cout << "x_ = " << x_ << endl;
	cout << "P_ = " << P_ << endl;



}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix, P_.
  */

	// 1. generate augmented sigma points; output: Xsig_aug
	MatrixXd Xsig_aug = GenerateSigmaPoints();

	// 2. predict augmented sigma points; output: set Xsig_pred_
	PredictSigmaPoints(Xsig_aug, delta_t);

	// 3. get state mean and covariance, set x_ and P_
	PredictNewState();

}

void UKF::PredictNewState() {
	//predict state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		x_ += weights_(i) * Xsig_pred_.col(i);
	}
	//predict state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P_ += weights_(i) * x_diff * x_diff.transpose();
	}
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_aug,double delta_t) {
	
	// set Xsig_pred_;
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd col = Xsig_aug.col(i);
		double px = col(0);
		double py = col(1);
		double v = col(2);
		double yaw = col(3);
		double yawd = col(4);
		double std_a = col(5);
		double std_yawdd = col(6);

		VectorXd xk = VectorXd(n_x_);
		xk << px,
			py,
			v,
			yaw,
			yawd;

		VectorXd noise = VectorXd(n_x_);
		noise << 0.5 * delta_t * delta_t * cos(yaw) * std_a,
			0.5 * delta_t * delta_t * sin(yaw) * std_a,
			delta_t * std_a,
			0.5 * delta_t * delta_t * std_yawdd,
			delta_t * std_yawdd;

		VectorXd transform = VectorXd(n_x_);
		if (fabs(yawd) <= 0.001) { // effectively zero
			transform << v * cos(yaw) * delta_t,
				v * sin(yaw) * delta_t,
				0,
				yawd * delta_t,
				0;
		}
		else {
			transform << v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw)),
				v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)),
				0,
				yawd * delta_t,
				0;
		}

		VectorXd next_state = xk + transform + noise;

		Xsig_pred_.col(i) << next_state;
	}
		
}

MatrixXd UKF::GenerateSigmaPoints() {
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug.head(n_x_) = x_;
	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;
	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();
	//create augmented sigma points
	Xsig_aug.col(0) << x_aug; 
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) << x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) << x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

	return Xsig_aug;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

	// normal kalman filter

	VectorXd z = meas_package.raw_measurements_;

	VectorXd z_pred = H_laser_ * x_; // dimension: 2 by 1
	VectorXd y = z - z_pred; // dimension: 2 by 1

	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;

	// update NIS_laser_
	
	NIS_laser_ = y.transpose() * Si * y;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

	// 1. get predicted sigma in measurement space
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd col = Xsig_pred_.col(i);
		double px = col(0);
		double py = col(1);
		double v = col(2);
		double yaw = col(3);
		double yawd = col(4);

		// divide by zero???
		if (abs(px) <= 0.0001 && abs(py) <= 0.0001) {
			px = 0.1;
			py = 0.1;
		}

		VectorXd z = VectorXd(n_z);
		z << sqrt(px*px + py*py),
			atan2(py, px),
			(px * cos(yaw)*v + py*sin(yaw)*v) / sqrt(px*px + py*py);

		Zsig.col(i) << z;
	}
	//calculate mean predicted measurement
	z_pred.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred += weights_(i) * Zsig.col(i);
	}
	//calculate measurement covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	S += R;

	// update measurement, update x_ and P_
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0],
	     meas_package.raw_measurements_[1],
	     meas_package.raw_measurements_[2];

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	Tc.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc += weights_(i) * x_diff * z_diff.transpose();

	}
	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	//update state mean and covariance matrix
	x_ += K * (z - z_pred);
	P_ = P_ - K * S * K.transpose();

	// update NIS_radar_
	VectorXd zdiff = z - z_pred;
	NIS_radar_ = zdiff.transpose() * (S.inverse()) * zdiff;

}
