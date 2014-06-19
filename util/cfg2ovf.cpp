/// Program to convert vampire cfg files to rasmol format
///
/// ./cfg2rasmol 

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <functional>

using namespace std;
using namespace std::tr2::sys;

typedef vector<path> vec;             // store paths,
                               // so we can sort them later
struct atom_coord_line_t {
    double cx, cy, cz;
    int mat, cat;
    std::string atom_type;
    atom_coord_line_t(double cx0, double cy0, double cz0,
        int mat0, int cat0, string atom_type0) :
        cx(cx0), cy(cy0), cz(cz0), mat(mat0), cat(cat0),
        atom_type(atom_type0)
    {};
};
typedef vector<atom_coord_line_t> coordLine_t;

struct spin_t {
    double sx, sy, sz;
    spin_t(double s1, double s2, double s3) :
        sx(s1), sy(s2), sz(s3) {};
};
typedef vector<spin_t> spin_lines_t;

void getCfgFiles(string dir, vec& fileNames)
{
    fileNames.clear();
    try
    {
        path p(dir);   // p reads clearer than argv[1] in the following code
        //p = p.parent_path();

        if (!exists(p))  {   // does p actually exist?
            cout << p << " does not exist\n";
            exit(0);
        }
        if (is_directory(p))      // is p a directory?
        {

            directory_iterator end_itr;
            for (directory_iterator it(p); it != end_itr; ++it) {
                if (it->path().extension() == ".cfg")
                    fileNames.push_back(it->path());
            }

            if (fileNames.size()<1) {
                cout << "No .cfg files\n";
                exit(0);
            }
            sort(fileNames.begin(), fileNames.end(), greater<string>());             // sort, since directory iteration
            // is not ordered on some file systems
        }
        else
            cout << p << " exists, but is not a directory\n";

    }
    catch (const filesystem_error& ex)  {
        cout << ex.what() << '\n';
    }
}

void readCoordFile(string coordFileName, coordLine_t& coord_lines, 
    double& stepx, double& stepy, double& stepz)
{
    coord_lines.clear();
    // open coordinate file
    std::ifstream coord_file;
    coord_file.open(coordFileName);


    // check for open file
    if (!coord_file.is_open()){
        std::cerr << "Error! Coordinate file atoms-coords.cfg cannot be opened. Exiting" << std::endl;
        exit(1);
    }

    // read in file header
    std::string dummy;

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    // get number of atoms
    unsigned int n_atoms;
    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;	
    dummy.erase(dummy.begin(), dummy.begin()+17);
    n_atoms = atoi(dummy.c_str());

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    // number of subsidiary files
    getline(coord_file, dummy);

    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;

    unsigned int n_local_atoms;
    getline(coord_file, dummy);
    //std::cout << dummy << std::endl;	
    n_local_atoms = atoi(dummy.c_str());

    double cx, cy, cz;
    int mat, cat;
    std::string atom_type;

    coord_lines.reserve(n_local_atoms);

    double cx0, cy0, cz0;
    cx0 = cy0 = cz0 = -100; // coords always >=0
    stepx = stepy = stepz = 1000000;
    // finish reading master file coordinates
    for (int i = 0; i<n_local_atoms; i++){
        coord_file >> mat >> cat >> cx >> cy >> cz >> atom_type;
        double temp;
        temp = abs(cx-cx0);
        if ((temp>0.1) && temp<stepx) stepx = temp;
        cx0 = cx;
        temp = abs(cy-cy0);
        if ((temp>0.1) && temp<stepy) stepy = temp;
        cy0 = cy;
        temp = abs(cz-cz0);
        if ((temp>0.1) && temp<stepz) stepz = temp;
        cz0 = cz;

        coord_lines.emplace_back(cx, cy, cz, mat, cat, atom_type);
    }

    // close master file
    coord_file.close();
}

void write_ovf_file(string filename, coordLine_t coord_lines, spin_lines_t spin_lines,
    double stepx, double stepy, double stepz)
{
    // Open output file
    std::ofstream outfile;
    outfile.open(filename);

    // Output ovf file header
    outfile << "# OOMMF: irregular mesh v0.0\n";
    outfile << "## File: " << filename << "\n";
    outfile << "## Grid step: " 
        << stepx << "\t" << stepy << "\t" << stepz << "\n";
    outfile << "# x y z m_x m_y m_z\n";

    for (int i = 0; i<coord_lines.size(); ++i) {
        auto line = coord_lines[i];
        auto spin = spin_lines[i];
        int mat = line.mat;
        if (mat>4) mat = 4;
        double ratio = 1.0 - 0.2*mat;
        outfile << line.cx << "\t" << line.cy << "\t" << line.cz << "\t"
            << spin.sx*ratio <<"\t" << spin.sy*ratio  <<"\t"<< spin.sz*ratio <<"\n";
    }
    
    outfile.close();
}

void read_spin_file(string spin_file, spin_lines_t& spin_lines)
{
    // atoms-XXXXXXXX.cfg
    std::ifstream infile;
    infile.open(spin_file);

    // read number of mats. in this file
    string dummy;
    for (int i = 1; i<=13; ++i) // line 17 is mat number
        getline(infile, dummy);
    dummy.erase(dummy.begin(), dummy.begin()+20);
    int n_local_mats = atoi(dummy.c_str());
    for (int i = 1; i<=n_local_mats; ++i)
        getline(infile, dummy);
    
    // Ignore num of spin files
    // TODO: ****only applied for serial****
    getline(infile, dummy);
    getline(infile, dummy);
    getline(infile, dummy);
    // read number of atoms in this file
    getline(infile, dummy);
    int n_local_atoms = atoi(dummy.c_str());

    spin_lines.clear();
    spin_lines.reserve(n_local_atoms);

    double sx, sy, sz;
    for (int i = 0; i<n_local_atoms; i++){
        infile >> sx >> sy >> sz;
        spin_lines.emplace_back(sx, sy, sz);
    }

    infile.close();
}

int main(int argc, char *argv[], char *envp[])
{
    string dir = current_path<string>();
    cout << "Working directory: " << dir <<"\n";
    vec filenames;

    getCfgFiles(dir, filenames);

    coordLine_t coord_lines;
    double stepx, stepy, stepz;
    readCoordFile(filenames[0], coord_lines, stepx, stepy, stepz);

	// now read spin configuration files
	for(int file=1; file<filenames.size(); file++){
		// atoms-XXXXXXXX.cfg
        spin_lines_t spin_lines;
		read_spin_file(filenames[file],spin_lines);

        if (spin_lines.size()!=coord_lines.size()) continue; // atom number do not match
        
        auto ovffile = filenames[file].replace_extension("ovf");
        
        write_ovf_file(ovffile, coord_lines, spin_lines, stepx, stepy, stepz);

        cout << "Writing to ... " << ovffile << "\n";

	}

    // finished
	return 0;

}
