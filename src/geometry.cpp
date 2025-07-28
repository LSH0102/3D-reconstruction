#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
//#include "Eigen/EigenSolver.h"
#include "cube_array.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <cmath>
#include <nlopt.hpp>

struct EdgeKey
{
	EdgeKey(unsigned i0, unsigned i1) :i0_(i0), i1_(i1) {}

	bool operator==(const EdgeKey& _rhs)const
	{
		return i0_ == _rhs.i0_ && i1_ == _rhs.i1_;
	}

	unsigned i0_, i1_;
};

struct EdgeHash
{
	std::size_t operator()(const EdgeKey& key) const
	{
		std::size_t seed = 0;
		seed ^= key.i0_ + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= key.i1_ + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return std::hash<std::size_t>()(seed);
	}
};

typedef std::unordered_map<EdgeKey, unsigned, EdgeHash> EMap;
typedef typename EMap::const_iterator  EMapIterator;
EMap edge2vertex;

extern const int edgeTable[256];
extern const int tritable[256][2][17];

int add_vertex(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, unsigned int i0, unsigned int i1, Eigen::MatrixXd& vertices, int& num_vertices, EMap& edge2vertex)
{
	EMapIterator it = edge2vertex.find(EdgeKey(i0, i1));
	if (it != edge2vertex.end())
	{
		return it->second;
	}


	double s0 = abs(values(i0));
	double s1 = abs(values(i1));
	double t = s0 / (s0 + s1);
	num_vertices++;
	if (num_vertices > vertices.rows())
	{
		vertices.conservativeResize(vertices.rows() + 10000, Eigen::NoChange);
	}
	vertices.row(num_vertices - 1) = (1 - t) * points.row(i0) + t * points.row(i1);
	edge2vertex[EdgeKey(i0, i1)] = num_vertices - 1;
	return num_vertices - 1;
}

void add_cube(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, const unsigned corner[], Eigen::MatrixXd& vertices, int& num_vertices, Eigen::MatrixXi& faces, int& num_faces, EMap& edge2vertex)
{
	unsigned char cubetype(0);
	int i;
	int samples[12];
	for (i = 0; i < 8; i++)
	{
		if (values[corner[i]] > 0.0)
		{
			cubetype |= (1 << i);
		}
	}
	if ((cubetype == 0) || (cubetype == 255))
	{
		return;
	}
	if (edgetable[cubetype] & 1)
		samples[0] = add_vertex(values, points, corner[0], corner[1], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 2)
		samples[1] = add_vertex(values, points, corner[1], corner[2], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 4)
		samples[2] = add_vertex(values, points, corner[3], corner[2], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 8)
		samples[3] = add_vertex(values, points, corner[0], corner[3], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 16)
		samples[4] = add_vertex(values, points, corner[4], corner[5], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 32)
		samples[5] = add_vertex(values, points, corner[5], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 64)
		samples[6] = add_vertex(values, points, corner[7], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 128)
		samples[7] = add_vertex(values, points, corner[4], corner[7], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 256)
		samples[8] = add_vertex(values, points, corner[0], corner[4], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 512)
		samples[9] = add_vertex(values, points, corner[1], corner[5], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 1024)
		samples[10] = add_vertex(values, points, corner[2], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 2048)
		samples[11] = add_vertex(values, points, corner[3], corner[7], vertices, num_vertices, edge2vertex);

	for (i = 0; tritable[cubetype][0][i] != -1; i += 3)
	{
		num_faces++;
		if (num_faces > faces.rows())
		{

			faces.conservativeResize(faces.rows() + 10000, Eigen::NoChange);
		}
		faces.row(num_faces - 1) <<
			samples[tritable[cubetype][0][i]],
			samples[tritable[cubetype][0][i + 1]],
			samples[tritable[cubetype][0][i + 2]];

	}
}

void Marching_Cubes(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, const int x_res, const int y_res, const int z_res, Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
	assert(points.cols() == 3);
	if (x_res < 2 || y_res < 2 || z_res < 2)
	{
		return;
	}
	int i;
	faces.resize(10000, 3);
	int n_faces = 0;
	vertices.resize(10000, 3);
	int n_vertices = 0;
	unsigned n_cubes = (x_res - 1) * (y_res - 1) * (z_res - 1);
	assert(unsigned(points.rows()) == x_res * y_res * z_res);

	unsigned int offsets[8];
	offsets[0] = 0;
	offsets[1] = 1;
	offsets[2] = 1 + x_res;
	offsets[3] = x_res;
	offsets[4] = x_res * y_res;
	offsets[5] = 1 + x_res * y_res;
	offsets[6] = 1 + x_res + x_res * y_res;
	offsets[7] = x_res + x_res * y_res;

	unsigned j;
	for (j = 0; j < n_cubes; ++j)
	{
		unsigned corner[8];
		for (i = 0; i < 8; ++i)
		{
			unsigned int idx = j;
			unsigned int X = x_res - 1, Y = y_res - 1;
			unsigned int x = idx % X; idx /= X;
			unsigned int y = idx % Y; idx /= Y;
			unsigned int z = idx;
			idx = x + y * x_res + z * x_res * y_res;
			corner[i] = idx + offsets[i];
		}
		add_cube(values, points, corner, vertices, n_vertices, faces, n_faces, edge2vertex);
	}
	vertices.conservativeResize(n_vertices, Eigen::NoChange);
	faces.conservativeResize(n_faces, Eigen::NoChange);

}


void compute_grad(int nx, int ny, int nz, double h, Eigen::SparseMatrix<double>& G)
{
	Eigen::SparseMatrix<double> Dx((nx - 1)* ny* nz, nx* ny* nz);
	Eigen::SparseMatrix<double> Dy(nx * (ny - 1) * nz, nx * ny * nz);
	Eigen::SparseMatrix<double> Dz(nx * ny * (nz - 1), nx * ny * nz);
	typedef Eigen::Triplet<double> T;
	std::vector<T> tx, ty, tz,t;
	int i, j, k;
	
	for (i = 0; i < nx-1; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				t.push_back(T(i + j * (nx-1)+k * (nx - 1) * ny, i + j * (nx - 1) + k * (nx - 1) * ny, -1.0 / h));
				t.push_back(T(i + j * (nx-1)+k * (nx-1)* ny, i + 1 + j * (nx - 1) + k * (nx - 1) * ny, 1.0 / h));
			}
		}

	}
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny-1; j++)
		{
			for (k = 0; k < nz; k++)
			{
				
				t.push_back(T(i + j * nx + k * nx * (ny-1)+ (nx - 1) * ny * nz, i + j * nx + k * nx * (ny - 1) ,-1.0 / h));
				t.push_back(T(i + j * nx + k * nx * (ny - 1) + (nx - 1) * ny * nz, i + (j + 1) * nx + k * nx * (ny - 1), 1.0 / h));
				
			}
		}

	}
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz-1; k++)
			{
				t.push_back(T(i + j * nx + k * nx * ny+ (nx - 1) * ny * nz+ nx * (ny - 1) * nz, i + j * nx + k * nx * ny, -1.0 / h));
				t.push_back(T(i + j * nx + k * nx * ny + (nx - 1) * ny * nz + nx * (ny - 1) * nz, i + j * nx + (k+1)*nx * ny, 1.0 / h));
				
			}
		}
		
	}
	
	G.setFromTriplets(t.begin(), t.end());
}

void poisson_recon(Eigen::MatrixXd &P, Eigen::MatrixXd &N, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	int n = P.rows();
	int nx, ny, nz;
	int i, j, k;
	double diam = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	const double pad = 8;
	double h = diam / (double)(30 + 2 * pad);
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + 2.0 * pad * h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + 2.0 * pad * h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + 2.0 * pad * h) / h, 3.);

	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				int ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}

	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx * ny * nz);
	Eigen::SparseMatrix<double> G((nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1), nx * ny * nz);
	compute_grad(nx, ny, nz, h, G);


	Eigen::SparseMatrix<double> Wx(P.rows(), (nx - 1) * ny * nz);
	Eigen::SparseMatrix<double> Wy(P.rows(), nx * (ny - 1) * nz);
	Eigen::SparseMatrix<double> Wz(P.rows(), nx * ny * (nz - 1));
	typedef Eigen::Triplet<double> T;
	std::vector<T>  tx, ty, tz, t;
	double alpha;
	int idx, idy, idz;
	for (i = 0; i < P.rows(); i++)
	{
		idx = (int)((P(i, 0) - corner(0) - h / 2) / h);
		idy = (int)((P(i, 1) - corner(1) - h / 2) / h);
		idz = (int)((P(i, 2) - corner(2) - h / 2) / h);
		alpha = (P(i, 0) - corner(0) - idx * h) / h;
		tx.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		tx.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, alpha));
		tx.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, alpha));
		tx.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, alpha));
		tx.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));

		alpha = (P(i, 1) - corner(1) - idy * h) / h;
		ty.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		ty.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, alpha));
		ty.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, alpha));
		ty.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, alpha));

		alpha = (P(i, 2) - corner(2) - idz * h) / h;

		tz.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));
	}
	Wx.setFromTriplets(tx.begin(), tx.end());
	Wy.setFromTriplets(ty.begin(), ty.end());
	Wz.setFromTriplets(tz.begin(), tz.end());
	Eigen::VectorXd vx = Wx.transpose() * N.col(0);
	Eigen::VectorXd vy = Wy.transpose() * N.col(1);
	Eigen::VectorXd vz = Wz.transpose() * N.col(2);
	Eigen::VectorXd v(vx.rows() + vy.rows() + vz.rows());
	v << vx,
		vy,
		vz;
	Eigen::SparseMatrix<double > L = G.transpose() * G;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver(L);
	v = G.transpose() * v;
	Eigen::VectorXd r = solver.solve(v);
	Eigen::VectorXd o = Eigen::VectorXd::Ones(P.rows());
	Eigen::VectorXd p = Eigen::VectorXd::Ones(r.rows());
	Eigen::SparseMatrix<double > W(P.rows(), r.rows());
	Eigen::VectorXd u(P.rows());
	double xd, yd, zd;
	double c00, c01, c10, c11, c0, c1, c;
	for (i = 0; i < P.rows(); i++)
	{
		idx = (int)((P(i, 0) - corner(0)) / h);
		idy = (int)((P(i, 1) - corner(1)) / h);
		idz = (int)((P(i, 2) - corner(2)) / h);
		xd = (P(i, 0) - corner(0) - idx * h) / h;
		yd = (P(i, 1) - corner(1) - idy * h) / h;
		zd = (P(i, 2) - corner(2) - idz * h) / h;

		c00 = r(idx + idy * nx + idz * nx * ny) * (1 - xd) + r(idx + 1 + idy * nx + idz * nx * ny) * xd;
		c01 = r(idx + idy * nx + (idz + 1) * nx * ny) * (1 - xd) + r(idx + 1 + idy * nx + (idz + 1) * nx * ny) * xd;
		c10 = r(idx + (idy + 1) * nx + idz * nx * ny) * (1 - xd) + r(idx + 1 + (idy + 1) * nx + idz * nx * ny) * xd;
		c11 = r(idx + (idy + 1) * nx + (idz + 1) * nx * ny) * (1 - xd) + r(idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny) * xd;
		c0 = c00 * (1 - yd) + c10 * yd;
		c1 = c01 * (1 - yd) + c11 * yd;
		c = c0 * (1 - zd) + c1 * zd;
		u(i) = c;
	}
	double sigma = o.dot(u) / n;

	g = r - sigma * p;
	Marching_Cubes(g, x, nx, ny, nz, V, F);
}

//code of modified Gauss reconstruction

void compute_w(Eigen::VectorXd &w, Eigen::MatrixXd Q,Eigen::MatrixXd P, double wmax,double wmin,int kw)
{
	int i, j, k;
	double x;
	double sum = 0.0;
	Eigen::MatrixXd R;
	Eigen::VectorXd v;
	for (i = 0; i < Q.rows(); i++)
	{
		R = P.rowwise() - Q.row(i);
		v=R.rowwise().squaredNorm();
		sum = 0.0;
		for (j = 0; j < kw; j++)
		{
			for (k = j; k < P.rows() - 1; k++)
			{
				if (v(k) > v(k + 1))
				{
					x = v(k + 1);
					v(k + 1) = v(k);
					v(k) = x;
				}
			}
			sum += v(j);
		}
		sum /= kw;
		
		sum = sqrt(sum);
		if (sum > wmax)
		{
			w(i) = wmax;
		}
		else if (sum < wmin)
		{
			w(i) = wmin;
		}
		else
		{
			w(i) = sum;
		}
	}
}


Eigen::MatrixXd compute_app(Eigen::MatrixXd P, Eigen::VectorXd w)
{
	int n = P.rows();
	int i, j, k;
	double r;
	Eigen::MatrixXd A(n, 3 * n);
	Eigen::MatrixXd phi(n, 3);
	Eigen::VectorXd v(3 * n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			r = (P.row(i) - P.row(j)).norm();
			if (r < w(i))
			{
				phi.row(j) = (P.row(j) - P.row(i)) / (4 * M_PI * w(i) * w(i) * w(i));
			}
			else
			{
				
				phi.row(j)= (P.row(j) - P.row(i)) / (4 * M_PI * r * r * r);
			}
		}
		//now we need to flatten the matrix phi to get a row vector A_{i}
		v << phi.col(0),
			phi.col(1),
			phi.col(2);
		A.row(i) = v.transpose();
	}
	return A;
}

Eigen::MatrixXd compute_aqp(Eigen::MatrixXd Q,Eigen::MatrixXd P, Eigen::VectorXd w)
{
	int n = Q.rows();
	int i, j, k;
	double r;
	Eigen::MatrixXd A(n, 3 * P.rows());
	Eigen::MatrixXd phi(P.rows(), 3);
	Eigen::VectorXd v(3 * P.rows());
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < P.rows(); j++)
		{
			r = (Q.row(i) - P.row(j)).norm();
			if (r < w(i))
			{
				phi.row(j) = (P.row(j) - Q.row(i)) / (4 * M_PI * w(i) * w(i) * w(i));
			}
			else
			{

				phi.row(j) = (P.row(j) - Q.row(i)) / (4 * M_PI * r * r * r);
			}
		}
		//now we need to flatten the matrix phi to get a row vector A_{i}
		v << phi.col(0),
			phi.col(1),
			phi.col(2);
		A.row(i) = v.transpose();
	}
	return A;
}

void gauss_recon(Eigen::MatrixXd& P, Eigen::MatrixXd& V, Eigen::MatrixXi& F,double wmin,double wmax,int kw)
{
	
	int n = P.rows();
	int nx, ny, nz;
	int i, j, k;
	double diam = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	const double pad = 8;
	double h = diam / (double)(30 + 2 * pad);
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + 2.0 * pad * h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + 2.0 * pad * h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + 2.0 * pad * h) / h, 3.);

	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				int ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}
	
	Eigen::VectorXd w(P.rows());
	compute_w(w, P, P,wmax, wmin, kw);
	Eigen::MatrixXd A(n, 3 * n);
	A = compute_app(P, w);
	double alpha = 1.05;
	Eigen::MatrixXd B = A * A.transpose();
	for (i = 0; i < n; i++)
	{
		B(i, i) += (alpha - 1.0) * A.row(i).norm()* A.row(i).norm();
	}
	Eigen::VectorXd b=0.5* Eigen::VectorXd::Ones(P.rows());
	Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(B);
	Eigen::VectorXd xi = solver.solve(b);

	

	Eigen::VectorXd mu = A.transpose() * xi;
	Eigen::VectorXd g(nx * ny * nz);
	
	compute_w(g, x, P, wmax, wmin, kw);
	
	Eigen::MatrixXd T = compute_aqp(x, P, g);
	
	Eigen::VectorXd f = T * mu;

	
	//set iso-value
	Eigen::VectorXd o = A * mu;
	
	Eigen::VectorXd e =o.mean() * Eigen::VectorXd::Ones(nx*ny*nz);
	f = f - e;

	Marching_Cubes(f, x, nx, ny, nz, V, F);
}


//code of variational reconstruction

void compute_M(Eigen::MatrixXd &M, Eigen::MatrixXd P)
{
	int n = P.rows();
	Eigen::MatrixXd M00(n, n);
	Eigen::MatrixXd M01(n, 3 * n);
	Eigen::MatrixXd M11(3*n, 3*n);
	int i, j, k,l;
	double r;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			r = (P.row(i) - P.row(j)).norm();
			M00(i, j) = r * r * r;
			M01(i, j) = -3 * r * (P(j, 0) - P(i, 0));
			M01(i, j+n) = -3 * r * (P(j, 1) - P(i, 1));
			M01(i, j+2*n) = -3 * r * (P(j, 2) - P(i, 2));
			for (k = 0; k < 3; k++)
			{
				for (l = 0; l < 3; l++)
				{
					if (r <1e-8)
					{
						M11( i +n* k, j +n* l) = 0.0;
					}
					else
					{
						
						if (l == k)
						{
							M11( i +n* k,  j +n* l) = -3 / r * (P(i, k) - P(j, k)) * (P(i, l) - P(j, l))-3 * r ;
						}
						else
						{
							M11(i + n * k, j + n * l) = -3 / r * (P(i, k) - P(j, k)) * (P(i, l) - P(j, l));
						}
					}
					
				}
			}
		}
	}
	M << M00, M01,
		M01.transpose(), M11;

}



void optim(Eigen::VectorXd& x, Eigen::VectorXd& grad, Eigen::VectorXd& m0, Eigen::VectorXd& v0, int& t, Eigen::MatrixXd H, int n, double lr)
{
	double beta1 = 0.9, beta2 = 0.999;
	double eps = 1e-8;
	t = t + 1;

	Eigen::VectorXd input(3 * n);
	Eigen::VectorXd sina_cosa_sinb_cosb(n * 4);
	int i, j, k;
	for (i = 0; i < n; ++i)
	{
		int ind = i * 4;
		sina_cosa_sinb_cosb(ind) = sin(x(2 * i));
		sina_cosa_sinb_cosb(ind + 1) = cos(x(2 * i));
		sina_cosa_sinb_cosb(ind + 2) = sin(x(2 * i + 1));
		sina_cosa_sinb_cosb(ind + 3) = cos(x(2 * i + 1));
	}
	for (i = 0; i < n; i++)
	{

		input(i) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 3);
		input(i + n) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 2);
		input(i + n * 2) = sina_cosa_sinb_cosb(4 * i + 1);
	}
	Eigen::VectorXd a2 = H * input;
	for (i = 0; i < n; i++)
	{
		grad(i * 2) = a2(i) * sina_cosa_sinb_cosb(4 * i + 1) * sina_cosa_sinb_cosb(4 * i + 3) + a2(i + n) * sina_cosa_sinb_cosb(4 * i + 1) * sina_cosa_sinb_cosb(4 * i + 2) - a2(i + n * 2) * sina_cosa_sinb_cosb(4 * i);
		grad(i * 2 + 1) = -a2(i) * sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 2) + a2(i + n) * sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 3);

	}
	m0 = beta1 * m0 + (1 - beta1) * grad;
	v0 = beta2 * v0 + (1 - beta2) * grad.cwiseProduct(grad);
	Eigen::VectorXd v1(2 * n);
	for (i = 0; i < 2 * n; i++)
	{
		v1(i) = (m0(i)) / (1 - pow(beta1, t)) / (sqrt(v0(i) / (1 - pow(beta2, t))) + eps);
	}
	x -= lr * v1;
	for (i = 0; i < 2 * n; i++)
	{
		if (x(i) > 2 * M_PI)
		{
			x(i) = M_PI;
		}
		if (x(i) < -2 * M_PI)
		{
			x(i) = -2 * M_PI;
		}
	}

}



void var_recon(Eigen::MatrixXd& P, Eigen::MatrixXd& V, Eigen::MatrixXi& F,double lambda)
{
	int n = P.rows();
	int nx, ny, nz;
	int i, j, k;
	double diam = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	const double pad = 8;
	double h = diam / (double)(30 + 2 * pad);
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + 2.0 * pad * h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + 2.0 * pad * h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + 2.0 * pad * h) / h, 3.);

	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				int ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}

	Eigen::MatrixXd M(4 * n, 4 * n);
	compute_M(M, P);


	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(4 * n, 4);
	for (i = 0; i < n; i++)
	{
		N(i, 0) = 1;
		for (j = 0; j < 3; j++)
		{
			N(i, j + 1) = P(i, j);
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 3; j++)
		{
			N(n + i + j * n, j + 1) = -1;
		}
	}

	Eigen::MatrixXd A(4 * n + 4, 4 * n + 4);
	A << M, N,
		N.transpose(), Eigen::MatrixXd::Zero(4, 4);

	Eigen::MatrixXd InvA = A.inverse();

	Eigen::MatrixXd J00 = InvA.block(0, 0, n, n);
	Eigen::MatrixXd J01 = InvA.block(0, n, n, 3 * n);
	Eigen::MatrixXd J11 = InvA.block(n, n, 3 * n, 3 * n);
	Eigen::MatrixXd B = (Eigen::MatrixXd::Identity(n, n) + lambda * J00).inverse();
	Eigen::MatrixXd H = J11 - lambda * J01.transpose() * B * J01;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
	Eigen::VectorXd D = es.eigenvalues();
	int min;
	double minv = D.minCoeff(&min);

	Eigen::VectorXd g = es.eigenvectors().col(min);





	Eigen::VectorXd init_normals(3 * n);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(4 * n);
	for (i = n; i < 4 * n; i++)
	{
		y(i) = g(i - n);
	}
	Eigen::Vector3d inim;
	for (i = 0; i < n; i++)
	{
		inim(0) = y(n + i);
		inim(1) = y(2 * n + i);
		inim(2) = y(3 * n + i);
		init_normals(3 * i) = y(n + i) / inim.norm();
		init_normals(3 * i + 1) = y(2 * n + i) / inim.norm();
		init_normals(3 * i + 2) = y(3 * n + i) / inim.norm();

	}

	Eigen::VectorXd solveval = Eigen::VectorXd::Zero(2 * n);
	Eigen::VectorXd grad(2 * n);
	//初始化角度参数
	for (i = 0; i < n; i++)
	{
		solveval(2 * i) = atan2(sqrt(init_normals(3 * i) * init_normals(3 * i) + init_normals(3 * i + 1) * init_normals(3 * i + 1)), init_normals(3 * i + 2));
		solveval(2 * i + 1) = atan2(init_normals(3 * i + 1), init_normals(3 * i));

	}
	Eigen::VectorXd m0 = Eigen::VectorXd::Zero(2 * n);
	Eigen::VectorXd v0 = Eigen::VectorXd::Zero(2 * n);
	int t = 0;
	for (i = 0; i < 5000; i++)
	{
		optim(solveval, grad, m0, v0, t, H, n, 0.1);
	}

	Eigen::VectorXd normals(3 * n);
	Eigen::VectorXd sina_cosa_sinb_cosb(n * 4);
	for (i = 0; i < n; i++)
	{
		int ind = i * 4;
		sina_cosa_sinb_cosb(ind) = sin(solveval(2 * i));
		sina_cosa_sinb_cosb(ind + 1) = cos(solveval(2 * i));
		sina_cosa_sinb_cosb(ind + 2) = sin(solveval(2 * i + 1));
		sina_cosa_sinb_cosb(ind + 3) = cos(solveval(2 * i + 1));
	}
	for (i = 0; i < n; i++)
	{

		normals(i) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 3);
		normals(i + n) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 2);
		normals(i + n * 2) = sina_cosa_sinb_cosb(4 * i + 1);
	}

	g = normals;
	

	Eigen::VectorXd s = -lambda * B * J01 * g;

	Eigen::VectorXd q(4 * n + 4);
	q << s,
		g,
		Eigen::VectorXd::Zero(4);

	q = InvA * q;






	Eigen::VectorXd a = q.block(0,0, n,1);
	Eigen::VectorXd b = q.block(n, 0, 3 * n, 1);
	
	Eigen::VectorXd c = q.block(4 * n+1, 0, 3, 1);
	double d = q(4*n);
	Eigen::VectorXd f = d*Eigen::VectorXd::Ones(nx * ny * nz);
	double val;
	double r;
	for (i = 0; i < nx * ny * nz; i++)
	{
		for (j = 0; j < P.rows(); j++)
		{
			r = (x.row(i) - P.row(j)).norm();
			f(i) += a(j) * r * r * r;
			f(i) += -b(j) * 3 * r *  (P(j, 0) - x.row(i)(0)) - b( j + n) * 3 * r *(P(j, 1) - x.row(i)(1)) - b(j + 2*n) * 3 * r *(P(j, 2) - x.row(i)(2));
			
		}
		f(i) += c.dot(x.row(i));

	}
	Marching_Cubes(f, x, nx, ny, nz, V, F);
}



void estimate_normals(Eigen::MatrixXd& P, Eigen::MatrixXd &N2,double lambda)
{
	int n = P.rows();
	int nx, ny, nz;
	int i, j, k;
	double diam = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	const double pad = 8;
	double h = diam / (double)(30 + 2 * pad);
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + 2.0 * pad * h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + 2.0 * pad * h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + 2.0 * pad * h) / h, 3.);

	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				int ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}

	Eigen::MatrixXd M(4 * n, 4 * n);
	compute_M(M, P);
	
	
	Eigen::MatrixXd N=Eigen::MatrixXd::Zero(4 * n, 4);
	for (i = 0; i < n; i++)
	{
		N(i, 0) = 1;
		for (j = 0; j <3; j++)
		{
			N(i, j + 1) = P(i, j);
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 3; j++)
		{
			N(n + i + j * n, j + 1) = -1;
		}
	}

	Eigen::MatrixXd A(4 * n + 4, 4 * n + 4);
	A << M, N,
		N.transpose(), Eigen::MatrixXd::Zero(4, 4);

	Eigen::MatrixXd InvA = A.inverse();

	Eigen::MatrixXd J00 = InvA.block(0, 0, n, n);
	Eigen::MatrixXd J01 = InvA.block(0, n, n, 3 * n);
	Eigen::MatrixXd J11 = InvA.block(n, n, 3 * n, 3 * n);
	Eigen::MatrixXd B = (Eigen::MatrixXd::Identity(n, n) + lambda * J00).inverse();
	Eigen::MatrixXd H = J11 - lambda * J01.transpose() * B * J01;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
	Eigen::VectorXd D = es.eigenvalues();
	int min;
	double minv = D.minCoeff(&min);

	Eigen::VectorXd g = es.eigenvectors().col(min);




	
	Eigen::VectorXd init_normals(3 * n);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(4 * n);
	for (i = n; i < 4 * n; i++)
	{
		y(i) = g(i - n);
	}
	Eigen::Vector3d inim;
	for (i = 0; i < n; i++)
	{
		inim(0) = y(n + i);
		inim(1) = y(2 * n + i);
		inim(2) = y(3 * n + i);
		init_normals(3 * i) = y(n + i) / inim.norm();
		init_normals(3 * i + 1) = y(2 * n + i) / inim.norm();
		init_normals(3 * i + 2) = y(3 * n + i) / inim.norm();

	}
	
	Eigen::VectorXd solveval= Eigen::VectorXd::Zero(2 * n);
	Eigen::VectorXd grad(2 * n);
	//初始化角度参数
	for (i = 0; i < n; i++)
	{
		solveval(2 * i) = atan2(sqrt(init_normals(3 * i) * init_normals(3 * i) + init_normals(3 * i + 1) * init_normals(3 * i + 1)), init_normals(3 * i + 2));
		solveval(2 * i + 1) = atan2(init_normals(3 * i + 1), init_normals(3 * i));

	}
	Eigen::VectorXd m0 = Eigen::VectorXd::Zero(2 * n);
	Eigen::VectorXd v0 = Eigen::VectorXd::Zero(2 * n);
	int t = 0;
	for (i = 0; i < 5000; i++)
	{
		optim(solveval, grad,m0,v0, t,H, n,0.1);
	}

	Eigen::VectorXd normals(3 * n);
	Eigen::VectorXd sina_cosa_sinb_cosb(n * 4);
	for (i = 0; i < n; i++)
	{
		int ind = i * 4;
		sina_cosa_sinb_cosb(ind) = sin(solveval(2 * i));
		sina_cosa_sinb_cosb(ind + 1) = cos(solveval(2 * i));
		sina_cosa_sinb_cosb(ind + 2) = sin(solveval(2 * i + 1));
		sina_cosa_sinb_cosb(ind + 3) = cos(solveval(2 * i + 1));
	}
	for (i = 0; i < n; i++)
	{

		normals(i) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 3);
		normals(i + n) = sina_cosa_sinb_cosb(4 * i) * sina_cosa_sinb_cosb(4 * i + 2);
		normals(i + n * 2) = sina_cosa_sinb_cosb(4 * i + 1);
	}

	g = normals;
	for (i = 0; i < n; i++)
	{
		N2(i, 0) = g( i);
		N2(i, 1) = g( i+n);
		N2(i, 2) = g( 2*n+i);
	}

}



