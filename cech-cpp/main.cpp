#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <regex>
#include <map>
#include <queue>
#include <limits>
#include <random>
#include <numeric>

// GUDHI includes for Cech Complex
#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

// Eigen for matrix operations and MDS
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// CGAL header for the exact d-dimensional kernel
#include <CGAL/Epeck_d.h>
using Kernel = CGAL::Epeck_d<CGAL::Dimension_tag<2>>; // We will embed in 2D
using Point = Kernel::Point_d;
using Points = std::vector<Point>;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;

enum class PermutationType { RAND, EDGE, DIST };

struct HiCData {
    int resolution;
    Eigen::MatrixXd frequency_matrix;
    double min_bin;
    double max_bin;
};

struct PersistenceInfo {
    int dimension;
    double birth;
    double death;
    std::vector<Vertex_handle> birth_simplex;
    std::vector<Vertex_handle> death_simplex;
};

class HiCAnalyzer {
private:
    std::string input_file;
    std::string output_path;
    std::string output_name;
    int resolution;
    bool enable_permutation_; // New member to control permutation
    std::mt19937 rng;

public:
    // Constructor updated to accept the permutation flag
    HiCAnalyzer(const std::string& input, const std::string& path, 
                const std::string& name, int res, bool enable_perm) 
        : input_file(input), output_path(path), output_name(name), 
          resolution(res), enable_permutation_(enable_perm) {
        std::random_device rd;
        rng.seed(rd());
    }

    HiCData readRawHiCData() {
        std::ifstream file(input_file);
        if (!file.is_open()) throw std::runtime_error("Cannot open input file: " + input_file);
        
        std::vector<std::vector<double>> all_data;
        std::string line;
        while (std::getline(file, line)) {
            line = std::regex_replace(line, std::regex("\\s+"), " ");
            std::istringstream iss(line);
            std::vector<double> row;
            double value;
            while (iss >> value) row.push_back(value);
            if (row.size() >= 3) all_data.push_back(row);
        }
        file.close();
        if (all_data.empty()) throw std::runtime_error("No valid data found in input file");

        double temp_min = std::numeric_limits<double>::max();
        double temp_max = std::numeric_limits<double>::lowest();
        for (const auto& row : all_data) {
            temp_min = std::min({temp_min, row[0], row[1]});
            temp_max = std::max({temp_max, row[0], row[1]});
        }

        int dim = static_cast<int>((temp_max - temp_min) / resolution + 1);
        std::cout << "Hi-C matrix size = " << dim << std::endl;
        Eigen::MatrixXd freq_matrix = Eigen::MatrixXd::Zero(dim, dim);
        for (const auto& row : all_data) {
            int i = static_cast<int>((row[0] - temp_min) / resolution);
            int j = static_cast<int>((row[1] - temp_min) / resolution);
            if (i >= 0 && i < dim && j >= 0 && j < dim) {
                freq_matrix(i, j) = row[2];
                freq_matrix(j, i) = row[2];
            }
        }

        HiCData result;
        result.resolution = resolution;
        result.frequency_matrix = freq_matrix;
        result.min_bin = temp_min;
        result.max_bin = temp_max;
        return result;
    }

    Eigen::MatrixXd normalizeMatrix(const Eigen::MatrixXd& tad_matrix) {
        double coeff = 100.0;

        // Compute log(x + 1) for all elements
        Eigen::MatrixXd log_matrix = (tad_matrix.array() + 1.0).log();

        // Compute log(max + 1)
        double log_max = std::log(tad_matrix.maxCoeff() + 1.0);

        // Apply the formula: log(max + 1) * 101 - log(x + 1) * 100
        Eigen::MatrixXd distance_matrix = log_max * (coeff + 1.0) - log_matrix.array() * coeff;

        return distance_matrix;
    }

    Points convertMatrixToPoints(const Eigen::MatrixXd& distance_matrix) {
        int n = distance_matrix.rows();
        if (n == 0) return {};
        Eigen::MatrixXd D2 = distance_matrix.array().square();
        Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n) - (1.0/n) * Eigen::MatrixXd::Ones(n, n);
        Eigen::MatrixXd B = -0.5 * J * D2 * J;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(B);
        if (eigen_solver.info() != Eigen::Success) {
            throw std::runtime_error("MDS: Eigendecomposition failed.");
        }
        Eigen::VectorXd eigenvalues = eigen_solver.eigenvalues();
        Eigen::MatrixXd eigenvectors = eigen_solver.eigenvectors();
        int dim_to_embed = 2;
        Eigen::VectorXd top_eigenvalues = eigenvalues.tail(dim_to_embed);
        Eigen::MatrixXd top_eigenvectors = eigenvectors.rightCols(dim_to_embed);
        for(int i = 0; i < top_eigenvalues.size(); ++i) {
            if (top_eigenvalues(i) < 0) top_eigenvalues(i) = 0;
        }
        Eigen::MatrixXd coords = top_eigenvectors * top_eigenvalues.cwiseSqrt().asDiagonal();
        Points points;
        points.reserve(n);
        for (int i = 0; i < n; ++i) {
            std::vector<double> p_coords = {coords(i, 1), coords(i, 0)}; 
            points.emplace_back(dim_to_embed, p_coords.begin(), p_coords.end());
        }
        return points;
    }

    Eigen::MatrixXd randomlyPermuteDistMat(const Eigen::MatrixXd& distance_matrix, PermutationType flag) {
        int n = distance_matrix.rows();
        Eigen::MatrixXd perm_mat = Eigen::MatrixXd::Zero(n, n);
        if (flag == PermutationType::RAND) {
            std::vector<double> values;
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) values.push_back(distance_matrix(i, j));
            std::shuffle(values.begin(), values.end(), rng);
            int k = 0;
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
                perm_mat(i, j) = perm_mat(j, i) = values[k++];
            }
        } else if (flag == PermutationType::EDGE) {
            std::vector<int> p(n);
            std::iota(p.begin(), p.end(), 0);
            std::shuffle(p.begin(), p.end(), rng);
            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
                perm_mat(i, j) = distance_matrix(p[i], p[j]);
            }
        } else if (flag == PermutationType::DIST) {
            for (int k = 1; k < n; ++k) {
                Eigen::VectorXd diag = distance_matrix.diagonal(k);
                double mean = diag.mean();
                double stddev = std::sqrt((diag.array() - mean).square().sum() / diag.size());
                if (stddev < 1e-10) stddev = 1e-10;
                std::normal_distribution<> dist(mean, stddev);
                for (int i = 0; i < n - k; ++i) {
                    double val = std::max(0.0, dist(rng));
                    perm_mat(i, i + k) = perm_mat(i + k, i) = val;
                }
            }
        }
        return perm_mat;
    }

    void writeDistanceMatrix(const Eigen::MatrixXd& dist_mat, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) { std::cerr << "Warning: Cannot write distance matrix to " << filename << std::endl; return; }
        file << dist_mat.format(Eigen::IOFormat(Eigen::FullPrecision, 0, "\t"));
        file.close();
        std::cout << "Wrote distance matrix to " << filename << std::endl;
    }

    void writeSkeleton(const Simplex_tree& st, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) { std::cerr << "Warning: Cannot write skeleton to " << filename << std::endl; return; }
        for (auto simplex : st.complex_simplex_range()) {
            if (st.dimension(simplex) > 0) {
                bool first = true;
                for (auto vertex : st.simplex_vertex_range(simplex)) {
                    if (!first) file << "\t";
                    file << vertex; first = false;
                }
                file << "\n";
            }
        }
        file.close();
        std::cout << "Wrote simplex skeleton to " << filename << std::endl;
    }

    void writePersistenceInfo(const std::vector<PersistenceInfo>& pers_info, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) { std::cerr << "Warning: Cannot write persistence info to " << filename << std::endl; return; }
        for (const auto& info : pers_info) {
            file << info.dimension << "\t" << info.birth << "\t" << info.death << "\t[";
            for (size_t i = 0; i < info.birth_simplex.size(); ++i) file << (i > 0 ? "," : "") << info.birth_simplex[i];
            file << "]\t[";
            for (size_t i = 0; i < info.death_simplex.size(); ++i) file << (i > 0 ? "," : "") << info.death_simplex[i];
            file << "]\n";
        }
        file.close();
        std::cout << "Wrote persistence pairs to " << filename << std::endl;
    }
    
    void performTDA(const Eigen::MatrixXd& distance_matrix, const std::string& file_prefix, double max) {
        std::cout << "--- Analyzing: " << file_prefix << " ---" << std::endl;
        
        Points points = convertMatrixToPoints(distance_matrix);
        if (points.empty()) {
            std::cerr << "Error: Point cloud is empty, cannot create Cech complex for " << file_prefix << std::endl;
            return;
        }
        
        Gudhi::cech_complex::Cech_complex<Kernel, Simplex_tree> cech_complex_from_points(points, max);

        Simplex_tree simplex_tree;
        cech_complex_from_points.create_complex(simplex_tree, 2);

        if (!simplex_tree.is_empty()) {
            std::cout << "Done creating simplex tree" << std::endl;
            
            // writeSkeleton(simplex_tree, file_prefix + "_skeleton.txt");
            
            Persistent_cohomology persistence(simplex_tree);
            persistence.init_coefficients(11); 
            persistence.compute_persistent_cohomology(-1.0);
            
            std::vector<PersistenceInfo> full_pers_info;
            auto persistence_pairs = persistence.get_persistent_pairs();
            for (const auto& pair : persistence_pairs) {
                PersistenceInfo info;
                info.dimension = simplex_tree.dimension(std::get<0>(pair));
                info.birth = simplex_tree.filtration(std::get<0>(pair)) / max;
                info.death = (std::get<1>(pair) != simplex_tree.null_simplex()) ? simplex_tree.filtration(std::get<1>(pair)) / max : std::numeric_limits<double>::infinity();
                for (auto v : simplex_tree.simplex_vertex_range(std::get<0>(pair))) info.birth_simplex.push_back(v);
                if (std::get<1>(pair) != simplex_tree.null_simplex()) {
                    for (auto v : simplex_tree.simplex_vertex_range(std::get<1>(pair))) info.death_simplex.push_back(v);
                }
                full_pers_info.push_back(info);
            }
            writePersistenceInfo(full_pers_info, file_prefix + "_persisdiagram.txt");
        } else {
            std::cerr << "Error: Failed to create Cech complex for " << file_prefix << std::endl;
        }
        std::cout << "--- Finished analysis for: " << file_prefix << " ---" << std::endl;
    }

    void processPermutations(const Eigen::MatrixXd& distance_matrix, double max) {
        std::cout << "\nProcessing random permutation (randperm)..." << std::endl;
        Eigen::MatrixXd rand_perm_mat = randomlyPermuteDistMat(distance_matrix, PermutationType::RAND);
        std::string rand_prefix = output_path + output_name + "_randperm";
        writeDistanceMatrix(rand_perm_mat, rand_prefix + "_distmat.txt");
        performTDA(rand_perm_mat, rand_prefix, max);

        std::cout << "\nProcessing edge permutation (edgeperm)..." << std::endl;
        Eigen::MatrixXd edge_perm_mat = randomlyPermuteDistMat(distance_matrix, PermutationType::EDGE);
        std::string edge_prefix = output_path + output_name + "_edgeperm";
        writeDistanceMatrix(edge_perm_mat, edge_prefix + "_distmat.txt");
        performTDA(edge_perm_mat, edge_prefix, max);

        std::cout << "\nProcessing distance permutation (distperm)..." << std::endl;
        Eigen::MatrixXd dist_perm_mat = randomlyPermuteDistMat(distance_matrix, PermutationType::DIST);
        std::string dist_prefix = output_path + output_name + "_distperm";
        writeDistanceMatrix(dist_perm_mat, dist_prefix + "_distmat.txt");
        performTDA(dist_perm_mat, dist_prefix, max);
    }

    // run() method updated to check the permutation flag
    void run() {
        auto start_time = std::chrono::high_resolution_clock::now();
        std::cout << "Resolution: " << resolution << std::endl;
        HiCData hic_data = readRawHiCData();
        std::cout << "Generate distance matrix..." << std::endl;
        Eigen::MatrixXd distance_matrix = normalizeMatrix(hic_data.frequency_matrix);

        double max = std::log(hic_data.frequency_matrix.maxCoeff() + 1.0) * 101;

        std::string original_prefix = output_path + output_name;
        writeDistanceMatrix(distance_matrix, original_prefix + "_distmat.txt");
        std::cout << "Calculate persistent homology for original data..." << std::endl;
        performTDA(distance_matrix, original_prefix, max);

        // Conditionally run permutations based on the flag
        if (enable_permutation_) {
            std::cout << "\nPermutation enabled. Generating and processing permutations..." << std::endl;
            processPermutations(distance_matrix, max);
        } else {
            std::cout << "\nPermutation disabled by default. Use --enable-permutation flag to run." << std::endl;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "\nTotal time (ms): " << duration.count() << std::endl;
    }
};

// Main function rewritten for robust argument parsing
int main(int argc, char* argv[]) {
    std::string input_file, output_name, output_path;
    int resolution = 0;
    bool enable_permutation = false; // Permutation is disabled by default

    const std::string usage = "Usage: " + std::string(argv[0]) + " -i <input_file> -o <output_name> -p <output_path> -r <resolution> [--enable-permutation]";
    
    int i = 1;
    while (i < argc) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            input_file = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            output_name = argv[++i];
        } else if (arg == "-p" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "-r" && i + 1 < argc) {
            try {
                resolution = std::stoi(argv[++i]);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid resolution value '" << argv[i] << "'. Must be an integer." << std::endl;
                return 1;
            }
        } else if (arg == "--enable-permutation") {
            enable_permutation = true;
        } else {
            std::cerr << "Error: Unknown or invalid argument '" << arg << "'" << std::endl;
            std::cerr << usage << std::endl;
            return 1;
        }
        i++;
    }
    
    if (input_file.empty() || output_name.empty() || output_path.empty() || resolution <= 0) {
        std::cerr << "Error: All required parameters (-i, -o, -p, -r) must be provided and valid." << std::endl;
        std::cerr << usage << std::endl;
        return 1;
    }
    
    if (!output_path.empty() && output_path.back() != '/') {
        output_path += '/';
    }
    
    try {
        HiCAnalyzer analyzer(input_file, output_path, output_name, resolution, enable_permutation);
        analyzer.run();
    } catch (const std::exception& e) {
        std::cerr << "An error occurred during analysis: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}