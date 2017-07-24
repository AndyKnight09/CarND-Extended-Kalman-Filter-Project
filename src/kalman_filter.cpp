#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Predict the state and covariance based on transition model
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Update the state and covariance based on new measurment
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// Update the state and covariance based on new measurment

	// Predict measurement using non-linear measurement model
	VectorXd zPred = CalculatePredictedMeasurement(x_);

	// Calculate innovation
	VectorXd y = z - zPred;
	
	// Wrap angle for y(1) between +/- pi
	if (y(1) < -M_PI) y[1] += 2 * M_PI;
	else if (y(1) > M_PI) y[1] -= 2 * M_PI;

	cout << "y = " << y << endl;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::CalculatePredictedMeasurement(const VectorXd& x_state)
{
	VectorXd zPred(3);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float r = sqrt(px*px + py*py);

	zPred << r, atan2(py, px), (px*vx + py*vy) / r;

	return zPred;
}
