#include <igl/boundary_facets.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/copyleft/cgal/extract_feature.h>
#include <igl/invert_diag.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>



//#include <igl/internal_angles.h>
//#include <igl/sharp_edges.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>


#include <CGAL/Triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <boost/function_output_iterator.hpp>

#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <numeric>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef CGAL::Triangulation_2<K>                            Triangulation;
typedef Triangulation::Point                                Point;

typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                      Point_2;
typedef Kernel::Point_3                      Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SMP = CGAL::Surface_mesh_parameterization;

//IGL functions
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd H;


void viewMeshIGL(const char* filename)
{
	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();
}

//I need to find a way to pass the mesh in by reference and perform operations on the mesh.
Eigen::MatrixXd getMeanCurvature()
{
	Eigen::MatrixXd HN;
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	HN = -Minv * (L*V);
	H = HN.rowwise().norm(); //up to sign
	return H;
}

void view2DTriangulation(const char* filename)
{
	std::ifstream in(filename);
	std::istream_iterator<Point> begin(in);
	std::istream_iterator<Point> end;

	Triangulation t;
	t.insert(begin, end);

	CGAL::draw(t);
}

struct halfedge2edge
{
	halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
		: m_mesh(m), m_edges(edges)
	{}
	void operator()(const halfedge_descriptor& h) const
	{
		m_edges.push_back(edge(h, m_mesh));
	}
	const Mesh& m_mesh;
	std::vector<edge_descriptor>& m_edges;
};

struct SamplerVariables
{
	float Rs;
	float Rf;
};

SamplerVariables CalibrateSampler( float a, float b, float c)
{
	float x1, x2, discriminant, realPart, imaginaryPart;
	discriminant = b * b - 4 * a*c;

	if (discriminant > 0) {
		x1 = (-b + sqrt(discriminant)) / (2 * a);
		x2 = (-b - sqrt(discriminant)) / (2 * a);
		std::cout << "Roots are real and different." << std::endl;
		std::cout << "x1 = " << x1 << std::endl;
		std::cout << "x2 = " << x2 << std::endl;
	}

	else if (discriminant == 0) {
		std::cout << "Roots are real and same." << std::endl;
		x1 = (-b + sqrt(discriminant)) / (2 * a);
		x2 = x1;
		std::cout << "x1 = x2 =" << x1 << std::endl;
	}

	else {
		throw("The sampler solution is not correct");
		realPart = -b / (2 * a);
		imaginaryPart = sqrt(-discriminant) / (2 * a);
		std::cout << "Roots are complex and different." << std::endl;
		std::cout << "x1 = " << realPart << "+" << imaginaryPart << "i" << std::endl;
		std::cout << "x2 = " << realPart << "-" << imaginaryPart << "i" << std::endl;

	}

	float Rf = std::max(x1, x2);
	if (Rf < 0) {
		throw("The sampler solution is not correct as the max is negative");
	}

	float Rs = 2 *Rf*Rf / sqrt(3);

	SamplerVariables result = { Rs, Rf };
	return result;
}

template <typename T>
void print(T const &input)
{
	for (int i = 0; i < input.size(); i++) {
		std::cout << input.at(i) << ' ';
	}
	std::cout << std::endl;
}

template <typename T>
bool vectorIsIn(T &a, T &b) {
	for (auto a_elt : a) {
		if (std::find(b.begin(), b.end(), a_elt) == b.end()) {
			return false;
		}
	}
	return true;
}

bool onlyOnePathLeft(std::vector<bool> paths) {
	int count = std::accumulate(paths.begin(), paths.end(), 0);
	if (count == 1) {
		return true;
	}
	else {
		return false;
	}
}

int whichPathLeft(std::vector<bool> path_bools) {
	for (int i = 0; i < path_bools.size(); i++) {
		if (path_bools[i] == true)
			return i;
	}
}

void checkProcessedOrBoundary(std::vector<int> &adjacent_face_list, std::vector<int> &bl_and_processed_vertices,
	std::vector<bool> &isBoundaryOrProcessed, Eigen::MatrixXd mesh_face_colour) 
{
	for (int i = 0; i < adjacent_face_list.size(); i++) {
		Eigen::VectorXi temp_vec;
		temp_vec = F.row(adjacent_face_list[i]);
		std::cout << temp_vec << std::endl;

		std::vector<int> unprocessed_adjacent_face_vertices(temp_vec.data(), temp_vec.data() + temp_vec.rows() * temp_vec.cols());

		print<std::vector<int>>(unprocessed_adjacent_face_vertices);

		//we need only two vertices to be present for it to be an edge and we know that there are only 3 vertice
		std::vector<int> comb_a = { unprocessed_adjacent_face_vertices[0], unprocessed_adjacent_face_vertices[1] };
		std::vector<int> comb_b = { unprocessed_adjacent_face_vertices[0], unprocessed_adjacent_face_vertices[2] };
		std::vector<int> comb_c = { unprocessed_adjacent_face_vertices[1], unprocessed_adjacent_face_vertices[2] };

		if (vectorIsIn(comb_a, bl_and_processed_vertices) || vectorIsIn(comb_b, bl_and_processed_vertices) || vectorIsIn(comb_c, bl_and_processed_vertices)) {
			isBoundaryOrProcessed[i] = 1;
			mesh_face_colour.row(adjacent_face_list[i]) << 255, 255, 255;
			std::cout << "The index that is definitely boundary is " << i << std::endl;

		}
		else {
			std::cout << "The index that is not boundary is " << i << std::endl;
			isBoundaryOrProcessed[i] = 0;
			mesh_face_colour.row(adjacent_face_list[i]) << 100, 0, 100;
		}
	}
}

void removeRow(Eigen::MatrixXd& matrix, int rowToRemove)
{
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows)
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

	matrix.conservativeResize(numRows, numCols);
}

std::vector<int> resamplingFeature(Eigen::MatrixXd &feature_edges, Eigen::MatrixXd & interpolated_edge_df, SamplerVariables &solution) {
	//Error diffusion for features
	//since the mesh has only boundary edges

	//Make a pseudo adjacency list
	std::vector <int> adjacent_edges = { 0 };
	int adjacency_edge_index = 0;
	for (int i = 0; i < feature_edges.rows(); i++) {

		int a = feature_edges(adjacency_edge_index, 0);
		int b = feature_edges(adjacency_edge_index, 1);

		for (int j = 0; j < feature_edges.rows(); j++) {
			if (a == feature_edges(j, 0) || a == feature_edges(j, 1) ||
				b == feature_edges(j, 0) || b == feature_edges(j, 1)) {

				if ((std::find(adjacent_edges.begin(), adjacent_edges.end(), j) != adjacent_edges.end()) || adjacency_edge_index == j) {
					//std::cout << "Are you always in???" << i << " " << j << std::endl;
				}
				else {
					adjacent_edges.push_back(j);
					adjacency_edge_index = j;
					break;
				}
			}
		}
	}

	std::cout << "adjacent_edges size" << adjacent_edges.size() << std::endl;
	std::cout << "number of feature_edges" << feature_edges.rows() << std::endl;
	print<std::vector<int>>(adjacent_edges);


	std::vector <int> no_of_samples_on_edge_vec(feature_edges.rows(), 0);

	//for the pseudo adjacency list
	int edge_index = adjacent_edges[0];
	float next_edge_error = 0;

	for (int i = 0; i < feature_edges.rows(); i++) {
		std::cout << "interpolated_edge_df for  " << edge_index << " is " << interpolated_edge_df(edge_index) << std::endl;
		std::cout << "next edge error for  " << edge_index << " is " << next_edge_error << std::endl;
		float edge_density = next_edge_error + interpolated_edge_df(edge_index);

		float num_samples_on_edge = edge_density * solution.Rf;
		std::cout << "num_samples_on_edge  " << num_samples_on_edge << std::endl;

		int num_samples_on_edge_round = round(num_samples_on_edge);
		std::cout << "num_samples_on_edge_round  " << num_samples_on_edge_round << std::endl;

		no_of_samples_on_edge_vec[edge_index] = num_samples_on_edge_round;

		next_edge_error = num_samples_on_edge - num_samples_on_edge_round;

		edge_index = adjacent_edges[i];
	}

	print<std::vector<int>>(no_of_samples_on_edge_vec);

	std::cout << "total number of samples on features is" << std::accumulate(no_of_samples_on_edge_vec.begin(), no_of_samples_on_edge_vec.end(), 0) << std::endl;
	return no_of_samples_on_edge_vec;
	//1. decide the backbones and the paths that can be traversed,
	//1a Look for vertices that are only present once in the whole feature array, they are possible starting points, choose any of them
	//1b Look for vertices that occur more than two times, i.e three times or more, those are corners, you can use this heuristic for Rs and Rf calculation as well.
	//Corners are to have a saample to themselves
	//for confluent corners, that have one edge already processed and passing errors to the other two edges, if you reach one of the processed edges again, add the error, 
	//process it, and get teleport the extra to the last feature. 
	//Find a way to chain the features that are close together and make them into a backbone
	//pick a vertex from 1a and start, use df/Rf to decide amount of samples on edge and place them there while recording them, and passing the error to the next feature
	//when you reach corners, pass them 
	//repeat process untill you are finished.
}

void surface_resampling(Eigen::MatrixXd &interpolated_face_ds, SamplerVariables &solution) {

	//set colour of faces
	Eigen::MatrixXd mesh_face_colour(F.rows(), 3);
	for (int i = 0; i < F.rows(); ++i) {
		mesh_face_colour.row(i) << 255, 255, 0;
	} // The idea behind this is to first set all the faces to green, then change all the faces processed to red or any other color

	std::vector <float> diffused_error(F.rows(), 0.0); // contains all the errors diffused to each face, originally set to zero
	std::vector <int> processed_faces = {};  //list of processed faces

	Eigen::MatrixXd bf, bf_into_F_indices, indices_of_vertex_across_bf;
	igl::boundary_facets(F, bf, bf_into_F_indices, indices_of_vertex_across_bf);

	std::cout << "did it work" << std::endl;
	std::cout << bf_into_F_indices(0) << std::endl;

	//Make a bf_into_F_indices_stdvec
	std::vector<int> bf_into_F_indices_stdvec(bf_into_F_indices.data(), bf_into_F_indices.data() + bf_into_F_indices.rows() * bf_into_F_indices.cols());
	std::cout << bf_into_F_indices_stdvec[0] << std::endl;

	std::vector <int> bl, bl_and_processed_vertices;
	igl::boundary_loop(F, bl);

	bl_and_processed_vertices = bl;
	std::cout << "boundary loop" << std::endl;
	/*print<std::vector<int>>(bl);*/

	int face_index = bf_into_F_indices(0); //select the first boundary face

	Eigen::MatrixXd TT, TTi, EdgeLengths;

	//get the triange adjancency
	igl::triangle_triangle_adjacency(F, TT, TTi);
	igl::edge_lengths(V, F, EdgeLengths);

	// Fix mis-match convention from igl
	Eigen::PermutationMatrix<3, 3> perm(3);
	perm.indices() = Eigen::Vector3i(1, 2, 0);
	TT = (TT*perm).eval();
	TTi = (TTi*perm).eval();
	for (int i = 0; i < TTi.rows(); i++) {
		for (int j = 0; j < TTi.cols(); j++) {
			TTi(i, j) = TTi(i, j) == -1 ? -1 : (int(TTi(i, j) + 3 - 1)) % 3;
		}
	}

	Eigen::VectorXd Area_values;
	igl::doublearea(V, F, Area_values);

	float total_area = (Area_values / 2).sum();
	float total_perimeter = EdgeLengths.sum();
	float compacity = total_area / total_perimeter;

	//Perform resampling
	for (int resampling_index = 0; resampling_index < F.rows(); ++resampling_index) { //F.rows()

		std::cout << "current loop = " << resampling_index << std::endl;
		float face_density = interpolated_face_ds(face_index);
		float samples_on_face = (face_density + diffused_error[face_index]) * solution.Rs;
		int samples_on_face_round = round(samples_on_face);

		float error_signed = samples_on_face - samples_on_face_round;

		/*std::cout << "samples on face is " << samples_on_face << std::endl;
		std::cout << "TT(face_index, 0) " << TT(face_index, 0) << std::endl;
		std::cout << "face_index " << face_index << std::endl;*/


		//get adjacent triangles
		std::vector<int> adjacent_faces = { int(TT(face_index, 0)), int(TT(face_index, 1)), int(TT(face_index, 2)) };
		std::vector<int> sorted_adjacent_faces = adjacent_faces;
		std::vector<int> unprocessed_adjacent_faces;

		//check if adjacent triangle is processed or not
		std::vector<int> v(3);
		std::vector<int>::iterator it;

		std::sort(sorted_adjacent_faces.begin(), sorted_adjacent_faces.end());
		std::sort(processed_faces.begin(), processed_faces.end());

		it = std::set_difference(sorted_adjacent_faces.begin(), sorted_adjacent_faces.end(), processed_faces.begin(), processed_faces.end(), v.begin());

		v.resize(it - v.begin());

		std::cout << "The difference has " << (v.size()) << " elements:\n";

		for (it = v.begin(); it != v.end(); ++it) {
			std::cout << ' ' << *it;
			if (*it != -1) {
				unprocessed_adjacent_faces.push_back(*it);
			}
		}

		if (unprocessed_adjacent_faces.size() < 1) {
			std::cout << "The current loop is " << resampling_index << " loop:\n";
			break;
		}

		std::cout << "The adjacent_faces are " << adjacent_faces[0] << " " << adjacent_faces[1] << " " << adjacent_faces[2] << std::endl;
		std::cout << "The sorted_adjacent_faces are " << sorted_adjacent_faces[0] << " " << sorted_adjacent_faces[1] << " " << sorted_adjacent_faces[2] << std::endl;
		std::cout << "The unprocessed face has size and value ";
		print<std::vector<int>>(unprocessed_adjacent_faces);

		//distribute error into the unprocessed ones
		std::vector<int> pos_of_unp_face_in_adjacency_face(unprocessed_adjacent_faces.size(), 0);
		for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
			int unprocessed_face_index = unprocessed_adjacent_faces[i];

			for (int j = 0; j < adjacent_faces.size(); j++) {
				if (unprocessed_face_index == adjacent_faces[j]) {
					pos_of_unp_face_in_adjacency_face[i] = j;
				}
			}
		}

		std::cout << "The unprocessed face pos in adj vec ";
		print<std::vector<int>>(pos_of_unp_face_in_adjacency_face);

		//get edge positions and edge lengths,
		std::vector<float> corresponding_edge_lengths(unprocessed_adjacent_faces.size(), 0);
		for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
			std::cout << i << std::endl;
			int pos = pos_of_unp_face_in_adjacency_face[i];
			int edge_position_index = TTi(face_index, pos);
			corresponding_edge_lengths[i] = EdgeLengths(unprocessed_adjacent_faces[i], edge_position_index);
		}

		std::cout << "The corresponding_edge_lengths ";
		print<std::vector<float>>(corresponding_edge_lengths);
		std::cout << "edge_lengths test" << EdgeLengths.row(face_index) << std::endl;


		//distribute the error

		/*Eigen::VectorXd  eigen_vector_lengths = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(corresponding_edge_lengths.data(), corresponding_edge_lengths.size());*/
		Eigen::VectorXd eigen_vector_lengths(corresponding_edge_lengths.size());
		for (int i = 0; i < corresponding_edge_lengths.size(); i++) {
			eigen_vector_lengths(i) = corresponding_edge_lengths[i];
		}

		float sum_of_edge_lengths = eigen_vector_lengths.sum();
		Eigen::VectorXd corresponding_error_diffusion = eigen_vector_lengths * error_signed / sum_of_edge_lengths;

		for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
			int unprocessed_face_index = unprocessed_adjacent_faces[i];
			diffused_error[unprocessed_face_index] += corresponding_error_diffusion[i];
		}

		std::cout << "The diffused errors" << diffused_error[unprocessed_adjacent_faces[0]] << std::endl;

		std::cout << "unprocessed_adjacent_faces.size() " << unprocessed_adjacent_faces.size() << std::endl;

		//move to the next face using bfs and compacity


		//Paths
		//Give priority
		if (unprocessed_adjacent_faces.size() == 1) {
			std::cout << "There is only one possible path so we choose that path" << std::endl;
			//flag f as processed
			processed_faces.push_back(face_index);
			//add new vertices to boundary and processed vertices
			Eigen::VectorXi flagged_face_vertices;
			flagged_face_vertices = F.row(face_index);
			for (int i = 0; i < flagged_face_vertices.size(); i++) {
				bl_and_processed_vertices.push_back(flagged_face_vertices(i));
			}
			/*std::cout << "boundary loop and processed vertice" << std::endl;
			print<std::vector<int>>(bl_and_processed_vertices);*/

			//colour it differently
			mesh_face_colour.row(face_index) << 0, 255, 150;

			face_index = unprocessed_adjacent_faces[0];
		}
		else {
			std::cout << "There is more than one possible path so we have to check heuristics" << std::endl;
			//Heuristic 1 check if it is edges are boundary, 
			std::vector<bool> isBoundary(unprocessed_adjacent_faces.size(), false);
			for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
				Eigen::VectorXi temp_vec;
				temp_vec = F.row(unprocessed_adjacent_faces[i]);
				std::cout << temp_vec << std::endl;

				std::vector<int> unprocessed_adjacent_face_vertices(temp_vec.data(), temp_vec.data() + temp_vec.rows() * temp_vec.cols());

				print<std::vector<int>>(unprocessed_adjacent_face_vertices);

				//we need only two vertices to be present for it to be an edge and we know that there are only 3 vertice
				std::vector<int> comb_a = { unprocessed_adjacent_face_vertices[0], unprocessed_adjacent_face_vertices[1] };
				std::vector<int> comb_b = { unprocessed_adjacent_face_vertices[0], unprocessed_adjacent_face_vertices[2] };
				std::vector<int> comb_c = { unprocessed_adjacent_face_vertices[1], unprocessed_adjacent_face_vertices[2] };

				if (vectorIsIn(comb_a, bl_and_processed_vertices) || vectorIsIn(comb_b, bl_and_processed_vertices) || vectorIsIn(comb_c, bl_and_processed_vertices)) {
					isBoundary[i] = 1;
					mesh_face_colour.row(unprocessed_adjacent_faces[i]) << 255, 255, 255;
					std::cout << "The index that is definitely boundary is " << i << std::endl;

				}
				else {
					std::cout << "The index that is not boundary is " << i << std::endl;
					isBoundary[i] = 0;
					mesh_face_colour.row(unprocessed_adjacent_faces[i]) << 100, 0, 100;
				}

			}

			std::cout << "Is any of the adjacent faces on the boundary " << std::endl;
			print<std::vector<bool>>(isBoundary);

			if (onlyOnePathLeft(isBoundary)) {
				std::cout << "I checked for only one remaining path and yes, there is only one " << std::endl;
				//flag f as processed
				processed_faces.push_back(face_index);
				//add new vertices to boundary and processed vertices
				Eigen::VectorXi flagged_face_vertices;
				flagged_face_vertices = F.row(face_index);
				for (int i = 0; i < flagged_face_vertices.size(); i++) {
					bl_and_processed_vertices.push_back(flagged_face_vertices(i));
				}
				/*std::cout << "boundary loop and processed vertice" << std::endl;
				print<std::vector<int>>(bl_and_processed_vertices);*/
				//colour it differently
				mesh_face_colour.row(face_index) << 0, 255, 150;


				//change face_index
				int PathLeft = whichPathLeft(isBoundary);
				std::cout << "The index that of the path left and which is chosen is " << PathLeft << std::endl;
				face_index = unprocessed_adjacent_faces[PathLeft];
				//std::cout << "Path chosen is " << unprocessed_adjacent_faces[PathLeft] << std::endl;

				continue;
			}


			//Heuristic 2 check the adjacents of the adjacents, if there is only one, then that is a boundary, if there are two of them, not boundary.
			// choose the one with only one.

			std::vector<bool> isAdjacentSharingBoundary(unprocessed_adjacent_faces.size(), false); //at - adjacent_triangle
			for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
				int at_face_index = unprocessed_adjacent_faces[i];

				std::vector<int> at_all_adjacent_faces = { int(TT(at_face_index, 0)), int(TT(at_face_index, 1)), int(TT(at_face_index, 2)) };
				std::vector<int> at_processed_adjacent_faces = {};

				//remove current face_index from list, -1, and any processed or boundary face 
				//If only one face remains, that is boundary
				for (int j = 0; j < at_all_adjacent_faces.size(); j++) {
					std::vector <int> temp_face_check = { at_all_adjacent_faces[j] };
					if (at_all_adjacent_faces[j] == face_index || at_all_adjacent_faces[j] == -1 ||
						vectorIsIn(temp_face_check, bf_into_F_indices_stdvec) || vectorIsIn(temp_face_check, processed_faces)) {
						std::cout << "I removed the face_index or -1 the value that I removed was " << at_all_adjacent_faces[j] << std::endl;
					}
					else {
						at_processed_adjacent_faces.push_back(at_all_adjacent_faces[j]);
					}
				}

				if (at_processed_adjacent_faces.size() == 1) {
					isAdjacentSharingBoundary[i] = true;
				}
				else {
					isAdjacentSharingBoundary[i] = false;
				}
			}
			std::cout << "heuristic 2 does the adjacent share a boundary? " << std::endl;
			print<std::vector<bool>>(isAdjacentSharingBoundary);


			if (onlyOnePathLeft(isAdjacentSharingBoundary)) {
				std::cout << "Heuristic 2 I checked for only one remaining path and yes, there is only one " << std::endl;
				//flag f as processed
				processed_faces.push_back(face_index);
				//add new vertices to boundary and processed vertices
				Eigen::VectorXi flagged_face_vertices;
				flagged_face_vertices = F.row(face_index);
				for (int i = 0; i < flagged_face_vertices.size(); i++) {
					bl_and_processed_vertices.push_back(flagged_face_vertices(i));
				}
				/*std::cout << "boundary loop and processed vertice" << std::endl;
				print<std::vector<int>>(bl_and_processed_vertices);*/
				//colour it differently
				mesh_face_colour.row(face_index) << 0, 255, 150;


				//change face_index
				int PathLeft = whichPathLeft(isAdjacentSharingBoundary);
				std::cout << "Heuristic 2 The index that of the path left and which is chosen is " << PathLeft << std::endl;
				face_index = unprocessed_adjacent_faces[PathLeft];
				//std::cout << "Path chosen is " << unprocessed_adjacent_faces[PathLeft] << std::endl;

				continue;
			}


			//if all three checks fail, then use the compacity stuff
			//quick test, besides that is how I would write the last test anyway
			//flag f as processed
			processed_faces.push_back(face_index);
			//add new vertices to boundary and processed vertices
			Eigen::VectorXi flagged_face_vertices;
			flagged_face_vertices = F.row(face_index);
			for (int i = 0; i < flagged_face_vertices.size(); i++) {
				bl_and_processed_vertices.push_back(flagged_face_vertices(i));
			}
			mesh_face_colour.row(face_index) << 0, 255, 0;

			face_index = unprocessed_adjacent_faces[0];
			std::cout << "Any path chosen " << std::endl;
			continue;
		}

		/*total_area = store_area_remaining;
		total_perimeter = store_perimeter_remaining;
		compacity = store_compacity;*/
	}

	//viewer
	igl::opengl::glfw::Viewer viewer_surface_sampling;
	viewer_surface_sampling.data().set_mesh(V, F);
	viewer_surface_sampling.data().set_colors(mesh_face_colour);


	//add edges

	viewer_surface_sampling.launch();
}


std::string parametrize(const char* filename)
{
	std::ifstream in(filename);
	if (!in) {
		std::cerr << "Problem loading the input data" << std::endl;
	}

	SurfaceMesh sm;
	in >> sm;
	// A halfedge on the border
	halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
	// The 2D points of the uv parametrisation will be written into this map
	typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
	UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;

	typedef SMP::Two_vertices_parameterizer_3<SurfaceMesh>  Border_parameterizer;
	typedef SMP::LSCM_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;


	SMP::Error_code err = SMP::parameterize(sm, Parameterizer(), bhd, uv_map);
	if (err != SMP::OK) {
		std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
	}

	std::string parametrized_mesh_name = "data/result.off";
	std::ofstream out(parametrized_mesh_name);
	SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

	return parametrized_mesh_name;
}

void CVT() {}

void lifting() {}



int main(int argc, char* argv[])
{
	const char* filename = (argc > 1) ? argv[1] : "data/lilium_s.off";
	std::ifstream input(filename);

	igl::readOFF(filename, V, F); //V vertices, F faces

	Eigen::MatrixXd J;
	Eigen::MatrixXd K;

	//get ds and df
	Eigen::MatrixXd H = getMeanCurvature(); //ds and df

	Eigen::MatrixXd feature_edges; //
	double angle_threshold = igl::PI * 0.11;
	//I had to type cast the line 86 and 87 on the extract_feature.cpp to make
	igl::copyleft::cgal::extract_feature(V, F, angle_threshold, feature_edges);

	//Extract feature vertices and view
	std::vector<int> feature_vertice_index(feature_edges.data(), feature_edges.data() + feature_edges.rows() * feature_edges.cols());
	sort(feature_vertice_index.begin(), feature_vertice_index.end());
	feature_vertice_index.erase(unique(feature_vertice_index.begin(), feature_vertice_index.end()), feature_vertice_index.end());

	igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(V, F);

	Eigen::MatrixXd feature_vertices(feature_vertice_index.size(), 3);
	for (int i = 0; i < feature_vertice_index.size(); ++i) {
		feature_vertices.row(i) = V.row(feature_vertice_index[i]);
	}

	viewer.data().point_size = 10;
	viewer.data().add_points(feature_vertices, Eigen::RowVector3d(255, 0, 0));

	///////add edges////
	for (int i = 0; i < feature_edges.rows(); ++i)
		viewer.data().add_edges
		(
			V.row(feature_edges(i, 0)),
			V.row(feature_edges(i, 1)),
			Eigen::RowVector3d(255, 255, 0)
		);

	viewer.launch();

	//////////////////Error diffusion//////////////////////
	///////Calibrate the sampler
	//Specify number of vertices
	float new_vertice_num = 800;
	float num_of_corners = 0;

	//interpolated face density
	Eigen::MatrixXd interpolated_face_ds(F.rows(), 1);

	for (int i = 0; i < F.rows(); ++i) {
		interpolated_face_ds(i) = (H(F(i, 0)) + H(F(i, 1)) + H(F(i, 2))) / 3;
	}

	//Interpolated edge density
	Eigen::MatrixXd interpolated_edge_df(feature_edges.rows(), 1);
	for (int i = 0; i < feature_edges.rows(); ++i) {
		interpolated_edge_df(i) = (H(feature_edges(i, 0)) + H(feature_edges(i, 1))) / 2;
	}


	//Quadractic equation, solving equation 1 in paper
	//Make ds and df
	

	float P = interpolated_face_ds.sum();
	float q = interpolated_edge_df.sum();

	std::cout << "P(ds) size and sum is " << P << " "<< interpolated_face_ds.rows() << std::endl;
	std::cout << "Q(df) size and sum is " << q << " "<< interpolated_edge_df.rows() << std::endl;

	
	float a_quad = 2 * P / sqrt(3);
	float b_quad = q;
	float c_quad = num_of_corners - new_vertice_num; //C-V

	std::cout << "a is " << a_quad << std::endl;
	std::cout << "b  is " << b_quad << std::endl;
	std::cout << "c is " << c_quad << std::endl;

	SamplerVariables solution = CalibrateSampler(a_quad, b_quad, c_quad);

	std::cout << "Rs solution is " << solution.Rs << std::endl;
	std::cout << "Rf solution is " << solution.Rf << std::endl;

	//Error diffusion for features
	 std::vector <int> no_of_samples_on_edge_vec = resamplingFeature(feature_edges, interpolated_edge_df, solution);
	 print<std::vector<int>>(no_of_samples_on_edge_vec);
	
	//Error diffusion for surfaces in progress
	 surface_resampling(interpolated_face_ds, solution);

	

	//Perform Parametrization

	Mesh mesh;
	if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
		std::cerr << "Not a valid input file." << std::endl;
		return 1;
	}

	std::string parametrized_mesh_name = parametrize(filename);
	std::cout << parametrized_mesh_name << std::endl;

	const char* Parametrized_filename = "data/result.off";
	std::ifstream in(Parametrized_filename);
	if (!in) {
		std::cerr << "Problem loading the input data" << std::endl;
	}



	

	//double target_edge_length = 0.1;
	//unsigned int nb_iter = 3;

	//std::cout << "Split border...";

	//std::vector<edge_descriptor> border;
	//PMP::border_halfedges(faces(mesh),
	//	mesh,
	//	boost::make_function_output_iterator(halfedge2edge(mesh, border)));
	//PMP::split_long_edges(border, target_edge_length, mesh);

	//std::cout << "done." << std::endl;

	//std::cout << "Start remeshing of " << filename
	//	<< " (" << num_faces(mesh) << " faces,";

	//std::vector<face_descriptor> seed, patch;
	//Mesh::Property_map<face_descriptor, int> selected
	//	= mesh.add_property_map<face_descriptor, int>("f:selected", 0).first;
	//seed.push_back(*(faces(mesh).first));
	//selected[seed.front()] = true;
	//CGAL::expand_face_selection(seed, mesh, 5, selected, std::back_inserter(patch));

	//std::cout << " and patch of size " << patch.size() << std::endl;
	//PMP::isotropic_remeshing(patch,
	//	target_edge_length,
	//	mesh,
	//	PMP::parameters::number_of_iterations(nb_iter)
	//	.face_patch_map(selected)
	//	.protect_constraints(true)//i.e. protect border, here
	//);

	//std::ofstream out("out.off");
	//out << mesh;
	//std::cout << "Remeshing done." << std::endl;

	return 0;
}






////plotting boundary extra stuff 
//Extract boundary faces and view.
	//std::vector<int> b_vertice_index(b.data(), b.data() + b.rows() * b.cols());
	//sort(b_vertice_index.begin(), b_vertice_index.end());
	//b_vertice_index.erase(unique(b_vertice_index.begin(), b_vertice_index.end()), b_vertice_index.end());

	//igl::opengl::glfw::Viewer viewer_boundary;

	//viewer_boundary.data().set_mesh(V, F);

	//Eigen::MatrixXd b_vertices(b_vertice_index.size(), 3);
	//for (int i = 0; i < b_vertice_index.size(); ++i) {
	//	b_vertices.row(i) = V.row(b_vertice_index[i]);
	//}

	///*viewer.data().point_size = 1;*/
	//viewer_boundary.data().add_points(b_vertices, Eigen::RowVector3d(0, 0, 0));

	//viewer_boundary.launch();


//Extra stuff on gettign ds, df and Rs
////////////////////Error diffusion//////////////////////
//	///////Calibrate the sampler
//	////Get Rs and Rf
//	//// i need to solve this with interpolation instead.
//	//Specify number of vertices
//float new_vertice_num = 10000;
//float num_of_corners = 0;
//
//
//
////Quadractic equation, solving equation 1 in paper
//
////Make ds and df
//Eigen::MatrixXd df(feature_vertice_index.size(), 1);
//Eigen::MatrixXd ds(H.rows() - feature_vertice_index.size(), 1);
//
//
//int j = 0; // where j is the index of the feature density function.
//int k = 0; // where k is the index of the surface density function
///*for (int i = 0; i < H.rows(); ++i) {
//	if (j < df.rows()) {
//		if (i == feature_vertice_index[j]) {
//			df.row(j) = H.row(i);
//			j++;
//			continue;
//		}
//	}
//	if ( k < ds.rows())
//	{
//		ds.row(k) = H.row(i);
//		k++;
//	}
//}*/
//
///*float P = ds.sum();
//float q = df.sum();*/
//
///*std::cout << "P(ds) size and sum is " << P << " "<< ds.rows() << std::endl;
//std::cout << "Q(df) size and sum is " << q << " "<< df.rows() << std::endl;*/
//
//
////float a_quad = 2 * P / sqrt(3);
////float b_quad = q;
////float c_quad = num_of_corners - new_vertice_num; //C-V
//
///*std::cout << "a is " << a_quad << std::endl;
//std::cout << "b  is " << b_quad << std::endl;
//std::cout << "c is " << c_quad << std::endl;*/
//
///*SamplerVariables solution = CalibrateSampler(a_quad, b_quad, c_quad);
//
//std::cout << "Rs solution is " << solution.Rs << std::endl;
//std::cout << "Rf solution is " << solution.Rf << std::endl;*/

/////for choosing path

//if (unprocessed_adjacent_faces.size() == 1) {
//	face_index = unprocessed_adjacent_faces[0];
//	store_area_remaining = total_area - Area_values(unprocessed_adjacent_faces[0]);
//	store_perimeter_remaining = total_perimeter - EdgeLengths.row(unprocessed_adjacent_faces[0]).sum();
//	store_compacity = store_area_remaining / store_perimeter_remaining;
//}
//else {
//	for (int i = 0; i < unprocessed_adjacent_faces.size(); i++) {
//		float temp_area_remaining = total_area - Area_values(unprocessed_adjacent_faces[i]);
//		float temp_perimeter_remaining = total_perimeter - EdgeLengths.row(unprocessed_adjacent_faces[i]).sum();
//		float calc_compatibility = temp_area_remaining / temp_perimeter_remaining;
//
//		//if two of its edges are boundary, and that makes it one,
//		/*std::vector <int> unprocessed_adjacent_face
//		 IsAdjacentTriangleOnBoundary(bl, unprocessed_adjacent_faces[i], )*/
//
//		 //if it is adjacent to a processed triangle and that makes it one
//
//		 //else choose with compatibility
//
//		if ()
//
//			if (calc_compatibility > store_compacity) {
//				store_compacity = calc_compatibility;
//				store_area_remaining = temp_area_remaining;
//				store_perimeter_remaining = temp_perimeter_remaining;
//				face_index = unprocessed_adjacent_faces[i];
//			}
//
//	}
//}
//total_area = store_area_remaining;
//total_perimeter = store_perimeter_remaining;
//compacity = store_compacity;
//
//std::cout << "face index" << face_index << std::endl;
//std::cout << "Total area" << total_area << std::endl;
//std::cout << "Total perimeter" << total_perimeter << std::endl;
//std::cout << "compacity" << compacity << std::endl;


//
//std::cout << "face index" << face_index << std::endl;
//std::cout << "Total area" << total_area << std::endl;
//std::cout << "Total perimeter" << total_perimeter << std::endl;
//std::cout << "compacity" << compacity << std::endl;
//
//float store_compacity = 0;
//float store_area_remaining = 0;
//float store_perimeter_remaining = 0;


//Eigen::VectorXd Area_values;
//igl::doublearea(V, F, Area_values);

//float total_area = (Area_values / 2).sum();
//float total_perimeter = EdgeLengths.sum();
//float compacity = total_area / total_perimeter;


	////Get feature edges using thresholding and boundary points
	//Eigen::MatrixXd E, uE;
	//Eigen::VectorXi EMAP;
	//std::vector<std::vector<Eigen::MatrixXd::Scalar>> uE2E;
	//igl::unique_edge_map(F, E, uE, EMAP, uE2E);

	// Constrain edges with a dihedral angle over 60°
	/*typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
	EIFMap eif = get(CGAL::edge_is_feature, mesh);
	PMP::detect_sharp_edges(mesh, 10, eif);
	int sharp_counter = 0;
	for (edge_descriptor e : edges(mesh))
		if (get(eif, e)) {
			std::cout << "CGAL e" << e << std::endl;
			std::cout << "next" << std::endl;
			++sharp_counter;
		}
		else {
			std::cout << " No sharp edge" << std::endl;
		}

	std::cout << sharp_counter << " sharp edges" << std::endl;*/