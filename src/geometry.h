#include "Eigen/Core"
#include <nlopt.hpp>

void poisson_recon(Eigen::MatrixXd& P, Eigen::MatrixXd& N, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void gauss_recon(Eigen::MatrixXd& P, Eigen::MatrixXd& V, Eigen::MatrixXi& F, double wmin, double wmax, int kw);
void var_recon(Eigen::MatrixXd& P, Eigen::MatrixXd& V, Eigen::MatrixXi& F, double lambda);
void estimate_normals(Eigen::MatrixXd& P, Eigen::MatrixXd& N2, double lambda);