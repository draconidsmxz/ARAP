#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/fit_rotations.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_rhs.h>
#include <igl/ARAPEnergyType.h>
#include <igl/columnize.h>
#include <igl/covariance_scatter_matrix.h>
#include <igl/speye.h>
#include <igl/repdiag.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

Eigen::MatrixXd V, U;
Eigen::MatrixXd R;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::SparseMatrix<Eigen::MatrixXd::Scalar> L;
Eigen::SparseMatrix<double> K;
double anim_t = 0.0;
double anim_t_dir = 0.03;
Eigen::RowVector3d mid;
igl::min_quad_with_fixed_data<double> solver_data;

//NewRotation

Eigen::SparseMatrix<double> CSM;
//End of NewRotation

Eigen::MatrixXd get_row(Eigen::MatrixXd M, int i) {
	Eigen::MatrixXd res(1, M.cols());
	for (int j = 0; j < M.cols();j++) {
		res(0, j) = M(i, j);
	}
	return res;
}

void optimizeRotation(Eigen::SparseMatrix<Eigen::MatrixXd::Scalar> &L) {
	using namespace std;

	int dim = V.cols();
	int n = V.rows();
	Eigen::MatrixXd S(n*dim, dim);

	for (int i = 0; i < n; i++) {

		Eigen::MatrixXd s(dim, dim);
		for (int p = 0; p < dim; p++) {
			for (int q = 0; q < dim; q++) {
				s(p, q) = 0;
			}
		}
		for (int j = 0; j < n; j++) {
			if (L.coeff(i, j) <= 0) { continue; }
			else {
				
				s += L.coeff(i, j) * ((get_row(V,i) - get_row(V,j)).transpose() * (get_row(V,i) - get_row(V,j)));
				cout << s << endl <<i<< endl;
			}
		}
		for (int p = 0; p < dim; p++) {
			for (int q = 0; q < dim; q++) {
				S(i * dim + p, q) = s(p, q);
			}
		}
	}
	S /= S.array().abs().maxCoeff();
	//cout << S << endl;
	if (R.rows() == 2)
	{
		igl::fit_rotations_planar(S, R);
	}
	else
	{
		igl::fit_rotations(S, true, R);
	}
	cout << R.size() << endl;
}



void optimizeRotation2(const Eigen::MatrixXd & bc) {
	using namespace Eigen;
	using namespace std;
	for (int bi = 0; bi < bc.rows(); bi++)
	{
		U.row(b(bi)) = bc.row(bi);
	}
	const auto & Udim = U.replicate(3, 1);
	MatrixXd S = CSM * Udim;
	S /= S.array().abs().maxCoeff();
	const int Rdim = 3;
	igl::fit_rotations(S, true, R);
	
}



void optimizePosition(const Eigen::MatrixXd & bc) {
	using namespace Eigen;
	using namespace std;
	VectorXd Rcol;
/*
	MatrixXd RR(3, 360);
	for (int i = 0; i < 360; i++) {
		for (int j = 0; j < 3; j++) {
			RR(j,i) = (i % 3==j);
		}
	}*/
	//cout << RR;

	/*
	for (int r = 0; r < 120; r++)
	{
		cout << R.block(0, 3 * r, 3, 3) << endl;
	}*/


	int Rdim = 3;
	int num_rots = K.cols() / Rdim / Rdim;
	int n = V.rows();
	igl::columnize(R, num_rots, 2, Rcol);
	//cout << K.size() << endl << R.size()<<endl<< Rcol.size() << endl;
	VectorXd Bcol = -K * Rcol;
	for (int c = 0; c < Rdim; c++)
	{
		VectorXd Uc, Bc, bcc, Beq;
		Bc = Bcol.block(c*n, 0, n, 1);
		if (bc.size() > 0)
		{
			bcc = bc.col(c);
		}
		igl::min_quad_with_fixed_solve(
			solver_data,
			Bc, bcc, Beq,
			Uc);
		U.col(c) = Uc;
	}

	return;
}


void MyArapPrecomputation()
{
	using namespace std;
	using namespace Eigen;
	typedef typename Eigen::MatrixXd::Scalar Scalar;
	typedef SparseMatrix<Scalar> SparseMatrixS;
	igl::cotmatrix(V, F, L);
	SparseMatrix<double> Q = (-L).eval();
	//cout << L << endl;

	igl::arap_rhs(V, F, 3, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, K);

	igl::min_quad_with_fixed_precompute(
		Q, b, SparseMatrix<double>(), true, solver_data);

	/*optimizeRotation(L);
	for (int r = 0; r < 120; r++)
	{
		cout << R.block(0, 3*r, 3, 3) << endl;
	}*/

	igl::covariance_scatter_matrix(V, F, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, CSM);
	/*
	SparseMatrix<double> G_sum;
	cout << G_sum << endl << endl;
	igl::speye(V.rows(), G_sum);
	SparseMatrix<double> G_sum_dim;
	cout << G_sum_dim << endl << endl;
	igl::repdiag(G_sum, 3, G_sum_dim);
	cout << G_sum_dim << endl << endl;
	CSM = (G_sum_dim * CSM).eval();
	*/
}





template <
typename Derivedbc,
typename DerivedU>
void solve(const Eigen::PlainObjectBase<Derivedbc> & bc, Eigen::PlainObjectBase<DerivedU> & U) 
{

	for (int i = 0; i <= 20; i++) {
		optimizeRotation2(bc);
		//std::cout << L;
		optimizePosition(bc);
	}
	/*



	for (int bi = 0; bi < bc.rows(); bi++)
	{
		U.row(b(bi)) = bc.row(bi);
		//std::cout << bc.row(bi)<<std::endl;
	}*/

}

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
	using namespace Eigen;
	using namespace std;
	MatrixXd bc(b.size(), V.cols());
	MatrixXd bc2(b.size(), V.cols());
	for (int i = 0; i < b.size(); i++)
	{
		bc.row(i) = V.row(b(i));
		bc2.row(i) = V.row(b(i));
		switch (S(b(i)))
		{
		case 0:
		{
			double theta = 0.5*sin(anim_t*2.)*sin(anim_t*2.)*igl::PI;
			bc(i, 0) = bc2(i, 0) * cos(theta) + bc2(i, 2)*sin(theta);
			bc(i, 2) = -bc2(i, 0) * sin(theta) + bc2(i, 2)*cos(theta);
			//bc(i, 1) +=  sin(i+anim_t)/2.;
			
			/*const double r = 1;
			bc(i, 0) += r*sin(0.5*anim_t*2.*igl::PI);
			bc(i, 2) -= r + r*cos(igl::PI + 0.5*anim_t*2.*igl::PI);
			//bc(i, 2) *= sin(0.5*i/10*anim_t*2.*igl::PI);
			*/break;
		}
		case 1:
		{
			const double r = 1;
			bc(i, 1) += r + r*cos(igl::PI + 0.15*anim_t*2.*igl::PI);
			bc(i, 2) -= r*sin(0.15*anim_t*2.*igl::PI);
			break;
		}
		case 2:
		{
			const double r =0;
			bc(i, 2) += r + r*cos(igl::PI + 0.35*anim_t*2.*igl::PI);
			bc(i, 0) += r*sin(0.35*anim_t*2.*igl::PI);
			break;
		}
		default:
			break;
		}
	}
	solve(bc, U);

	//igl::arap_solve(bc, arap_data, U);
	viewer.data().set_vertices(U);
	viewer.data().compute_normals();
	
	if (viewer.core.is_animating)
	{
		anim_t += anim_t_dir;
	}
	return false;
}


bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
	switch (key)
	{
	case ' ':
		std::cout << "is_animating = " << !viewer.core.is_animating << std::endl;
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}


int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;
	// Load a mesh in OFF format
	igl::readOFF("D:/software/asap/Meshes_ARAP_SA2007/bar3.off", V, F);
	U = V;
	igl::readDMAT("D:/software/asap/Meshes_ARAP_SA2007/bar3_selection.dmat", S);

	// vertices in selection
	igl::colon<int>(0, V.rows() - 1, b);
	b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
		[](int i)->bool {return S(i) >= 0; }) - b.data());
	// Centroid
	mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());

	// Precomputation

	for (int i = 0; i < b.size(); i++) {
		cout << b[i] << " ";
	}
	cout << endl;

	MatrixXd C(F.rows(), 3);
	RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
	RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
	for (int f = 0; f < F.rows(); f++)
	{
		if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
		{
			C.row(f) = purple;
		}
		else
		{
			C.row(f) = gold;
		}
	}
	
	MyArapPrecomputation();

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.data().set_colors(C);
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;
	viewer.launch();
}

/*
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
*/