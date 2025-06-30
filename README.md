# topologicalhic

This repository contains tools for applying topological data analysis (TDA) to Hi-C data. The codebase includes both Python and C++ implementations for analyzing chromosome structure through persistent homology, with optional metric correction capabilities.

## 📋 Features

### Core Analysis Methods
- **Čech Complex Analysis** (C++)
- **Rips Complex Analysis** (Python + C++ binding)
- **Alpha Complex Analysis** (Python)
- **Metric Correction**: HLWB algorithm 


## 🛠 Installation

### Dependencies

#### C++ Implementation
```bash
# Required libraries
sudo apt-get install libboost-all-dev libeigen3-dev libcgal-dev libtbb-dev

# GUDHI library
git clone https://github.com/GUDHI/gudhi-devel.git
```

#### Python Implementation
```bash
pip install numpy scipy scikit-learn gudhi networkx

# For matrix correction (compile the C++ extension)
pip install pybind11
python setup.py build_ext --inplace
pip install .
```

### Building

#### C++ Version
```bash
mkdir build && cd build
cmake ..
make
```

#### Python Extension
```bash
python setup.py build_ext --inplace
pip install .
```

## 🚀 Usage

### Cech Complex C++ Implementation

```bash
# Basic analysis
./your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount

# With permutation analysis (disabled by default for performance)
./your_file_name -i input_data.txt -p output_path/  -o output_name -r resolution_amount --enable-permutation
```

### Vietoris-Rips Python Implementation

```bash
# Standard analysis
python your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount

# With metric correction
python your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount --correct_matrix

# With permutation controls
python your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount --run_permutations

# Full analysis with all options
python your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount --correct_matrix --run_permutations
```

### Alpha Complex

```bash
python your_file_name -i input_data.txt -p output_path/ -o output_name -r resolution_amount
```

### Sample
```bash
!python hic.py -i WTC_PSC_combined_200000_iced_chr15.matrix -p wtc_results/WTC_PSC/ -o WTC_PSC_combined_200000_iced_chr15 -r 200000
```

## 📊 Input Format

Hi-C contact data should be in tab/space-separated format:
```
#sample
genomic_position_1    genomic_position_2   contact_frequency
0                     200000               42.5
0                     400000               18.2
20000                 400000               35.7
...
```

## 📈 Output Files

### Standard Outputs
- `*_distmat.txt`: Distance matrix derived from contact frequencies
- `*_persisdiagram.txt`: Persistence diagram with homology dimensions, birth/death times with persitence pairs
  
```bash
#sample output 
RUES_CM_combined_200000_iced_chr14_distmat.txt
RUES_CM_combined_200000_iced_chr14_persisdiagram.txt
```

### Loop Information Output
- `*_loop_size_analysis.txt`: Information about loops such as loop_size_bins, genomic_size.
- `*_loop_size_analysis_summary.txt`: Summary statistics for loops
  
### Permutation Outputs (for human chr13-22)
- `*_randperm_persisdiagram.txt`: Random permutation analysis
- `*_edgeperm_persisdiagram.txt`: Edge permutation analysis  
- `*_distperm_persisdiagram.txt`: Distance-preserving permutation analysis

### Other Outputs
- `*_skeleton.txt`:  nodes forming a simplex and birth time of simplex

## 🔧 Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-i, --input` | Input Hi-C contact matrix file | Required |
| `-o, --output` | Base name for output files | Required |
| `-p, --path` | Output directory path | Required |
| `-r, --resolution` | Hi-C data resolution (bp) | Required |
| `--enable-permutation` | Enable permutation analysis (C++ cech complex) | Disabled |
| `--correct_matrix` | Apply HLWB metric correction (Python) | Disabled |
| `--run_permutations` | Run permutation controls (Python - vietoris rips) | Disabled |

* For permutations in the Alpha Complex - comment out the specific lines

## 🧮 Algorithms

### HLWB Metric Correction
A two-stage algorithm for projecting distance matrices to valid metric spaces:
1. **Stage I**: Kernel projection with exponential mapping
2. **Stage II**: Iterative triangle inequality enforcement

### Topological Analysis
- **Čech Complex**
- **Rips Complex**
- **Alpha Complex**

## 📖 Citation

If you use this toolkit in your research, please cite the appropriate methods:
- GUDHI library for topological data analysis
- HLWB algorithm for metric correction
- Hi-C data processing methodologies

## 🐛 Troubleshooting

### Common Issues
1. **Missing GUDHI**: Ensure GUDHI is properly installed and paths are correct
2. **Memory errors**: Large matrices may require substantial RAM
3. **Compilation errors**: Check all dependencies are installed
4. **Matrix correction failures**: Ensure input matrices are symmetric

### Performance Tips
- Use C++ version for production runs
- Enable permutations only when needed for statistical validation
- Use metric correction selectively for final analyses

## 📄 License

This project is provided as-is for research purposes. Please ensure compliance with the licenses of all dependencies (GUDHI, CGAL, etc.).
