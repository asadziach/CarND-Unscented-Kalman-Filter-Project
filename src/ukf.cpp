#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
	std_a_ = 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 30;

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

	// Parameters above this line are scaffolding, do not modify

	/**
	 TODO:

	 Complete the initialization. See ukf.h for other member properties.

	 Hint: one or more values initialized above might be wildly off...
	 */
	time_us_ = 0;

	//set state dimension
	n_x_ = 5;

	//set augmented dimension
	n_aug_ = 7;

	// TODO Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.2;

	//TODO Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.2;

	//define spreading parameter
	lambda_ = 3 - n_aug_;

	//x_aug_ = VectorXd(n_aug_);

	//P_aug_ = MatrixXd(n_aug_, n_aug_);

	int sig_pt_count = 2 * n_aug_ + 1;

	//create sigma point matrix
	//Xsig_aug_ = MatrixXd(n_aug_, sig_pt_count);

	//create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, sig_pt_count);

	//create vector for weights
	weights_ = VectorXd(sig_pt_count);

	//set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i < sig_pt_count; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	is_initialized_ = false;
}

UKF::~UKF() {
}

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
	if (!is_initialized_) {
		/**
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 */

		//state covariance matrix P

		P_   <<   1,  0,  0,  0,  0,
		          0,  1,  0,  0,  0,
		          0,  0,  1,  0,  0,
		          0,  0,  0,  1,  0,
		          0,  0,  0,  0,  1;

		x_.fill(0.0);

		float px,py;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			 */
			// Range - Radial distance from origin
			float rho = meas_package.raw_measurements_[0];
			// Bearing - angle between rho and x
			float phi = meas_package.raw_measurements_[1];
			// Radial Velocity - change of rho - range rate
			float rho_dot = meas_package.raw_measurements_[2];

			px = rho * cos(phi);
			py = rho * sin(phi);
			/*  although radar does include velocity information, the radar velocity
			 * and the CTRV velocity are not the same.
			 * */


		} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			 */
			px = meas_package.raw_measurements_[0];
			py = meas_package.raw_measurements_[1];

		}
		if (fabs(px) < 0.0001) {
			px = 0.01;
			cout << "initial px too small, clamping" << endl;
		}

		if (fabs(py) < 0.0001) {
			py = 0.01;
			cout << "initial py too small, clamping" << endl;
		}

		x_ << px, py, 0, 0, 0;

		time_us_ = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;

	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	//compute the time elapsed between the current and previous measurements. dt - expressed in seconds
	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	int sig_pt_count = 2 * n_aug_ + 1;

	//TOOD move up
	VectorXd x_aug_ = VectorXd(n_aug_);
	MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
	MatrixXd Xsig_aug_ = MatrixXd(n_aug_, sig_pt_count);

	//create augmented mean state
	x_aug_.head(n_x_) = x_;
	x_aug_ (n_x_) = 0;
	x_aug_(n_x_ + 1) = 0;

	//create augmented covariance matrix
	P_aug_.setZero();
	P_aug_.topLeftCorner(n_x_, n_x_) = P_;
	P_aug_(n_x_, n_x_) = std_a_ * std_a_;
	P_aug_(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd A = P_aug_.llt().matrixL();

	//create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;

	float m_const = sqrt(lambda_ + n_aug_);
	for (int i = 0; i < n_aug_; i++) {
		int k = i + 1;
		Xsig_aug_.col(k) = x_aug_ + m_const * A.col(i);
		Xsig_aug_.col(k + n_aug_) = x_aug_ - m_const * A.col(i);
	}

	Prediction(dt);

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
	 vector, x_. Predict sigma points, the state, and the state covariance matrix.
	 */
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
}
