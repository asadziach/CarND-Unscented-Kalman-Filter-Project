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

	// initial covariance matrix. For the CTRV model, P is a 5x5
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1.35;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.35;

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

	if (!is_initialized_) {
		/**
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 */

		float px,py;
		float variance_px, variance_py;
		float variance_v = pow(6.5, 2);
		float varience_psi = pow(M_PI, 2);
		float varience_psi_dot = pow(0.2, 2);

		P_.fill(0.0);
		x_.fill(0.0);

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
			float std_radr_max = max(std_radr_, std_radphi_ * rho);
			float variacne_radar = std_radr_* std_radr_;

			variance_px = variance_py = variacne_radar;

		} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			 */
			px = meas_package.raw_measurements_[0];
			py = meas_package.raw_measurements_[1];

			variance_px = std_laspx_* std_laspx_;
			variance_py= std_laspy_ * std_laspy_;

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

		// Set diagonal values
		P_(0, 0) = variance_px;
		P_(1, 1) = variance_py;
		P_(2, 2) = variance_v;
		P_(3, 3) = varience_psi;
		P_(4, 4) = varience_psi_dot;

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


	Prediction(dt);

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/
	/**
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	 */

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
		// Radar updates
		UpdateRadar(meas_package);
	} else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
		// Laser updates
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	 Estimate the object's location. Modify the state vector, x_.
	 Predict sigma points, the state, and the state covariance matrix.
	 */

	int sig_pt_count = 2 * n_aug_ + 1;

	//TODO try moving up
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

	//predict sigma points
	for (int i = 0; i < sig_pt_count; i++) {
		//extract values for better readability
		double p_x = Xsig_aug_(0, i);
		double p_y = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double yaw = Xsig_aug_(3, i);
		double yawd = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		} else {
			px_p = p_x + v * delta_t * cos(yaw);
			py_p = p_y + v * delta_t * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	//Update sate

	//predict state mean
	x_.fill(0.0);
	for (int i = 0; i < sig_pt_count; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	//predict state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < sig_pt_count; i++) {  //iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		float angle = x_diff(3);

		//angle normalization

		while (angle > M_PI)
			angle -= 2. * M_PI;
		while (angle < -M_PI)
			angle += 2. * M_PI;

		x_diff(3) = angle;

		//x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	 Use lidar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.
	 Finally calculate the lidar NIS.
	 */
	MatrixXd H_ = MatrixXd(2, 5);

	//Laser H matrix
	H_ << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0;

	//Sensor data

	float px = meas_package.raw_measurements_(0);
	float py = meas_package.raw_measurements_(1);

	VectorXd z = VectorXd(2);
	z << px, py;

	MatrixXd R_laser = MatrixXd(2, 2);
	R_laser << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_laser;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

	//Calculate Normalized Innovations Squared for LASAR
	float nis = y.transpose() * S.inverse() * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	 Use radar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.

	 Finally calculate the radar NIS.
	 */

	  //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//define spreading parameter
	double lambda = 3 - n_aug_;

	//set vector for weights
	VectorXd weights = VectorXd(2 * n_aug_ + 1);
	double weight_0 = lambda / (lambda + n_aug_);
	weights(0) = weight_0;
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {
		double weight = 0.5 / (n_aug_ + lambda);
		weights(i) = weight;
	}

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;
		double v2 = sin(yaw) * v;

		// measurement model
		Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                         //r
		Zsig(1, i) = atan2(p_y, p_x);                                     //phi
		Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0, std_radrd_
			* std_radrd_;
	S = S + R;

	//Update x_ and P_

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//Sensor data

	float rho = meas_package.raw_measurements_(0);
	float phi = meas_package.raw_measurements_(1);
	float rhodot = meas_package.raw_measurements_(2);

	VectorXd z = VectorXd(n_z);
	z << rho, phi, rhodot;

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization

		while (z_diff(1) > M_PI)
			z_diff(1) -= 2. * M_PI;
		while (z_diff(1) < -M_PI)
			z_diff(1) += 2. * M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1) > M_PI)
		z_diff(1) -= 2. * M_PI;
	while (z_diff(1) < -M_PI)
		z_diff(1) += 2. * M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

	//Calculate Normalized Innovations Squared for RADAR
	float nis = z_diff.transpose() * S.inverse() * z_diff;
}

