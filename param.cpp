#include "param.h"
#include <iostream>
#include <fstream>
#include <sstream>

template<typename T>
T getToken(const std::string& line)
{
	T	value;
    std::stringstream    ss(line);
    ss >> value;

	return value;
}

void    param::load_param(const std::string  filename)
{
	std::ifstream	infile(filename.c_str(),std::ifstream::in);	
	if (infile.good())
	{
		std::string line;

		while(getline(infile, line)) {

			std::string tmp, value;
			std::stringstream    ss(line);
		    ss >> tmp >> value;

			if (tmp.compare("part_input_file") == 0) {
				part_input_file = getToken<std::string>(value); 
			} else if (tmp.compare("dt") == 0) {
				dt = getToken<cl_real>(value); 
			} else if (tmp.compare("time_end") == 0) {
				time_end = getToken<cl_real>(value);
			} else if (tmp.compare("h") == 0) {
				h = getToken<cl_real>(value);
			} else if (tmp.compare("XCells") == 0) {
				XCells = getToken<cl_uint>(value);
			} else if (tmp.compare("YCells") == 0) {
				YCells = getToken<cl_uint>(value);
			} else if (tmp.compare("ZCells") == 0) {
				ZCells = getToken<cl_uint>(value);
			} else if (tmp.compare("part_out_freq") == 0) {
				part_out_freq = getToken<cl_uint>(value);
			} else if (tmp.compare("part_out_name_base") == 0) {
				part_out_name_base = getToken<std::string>(value);
			} else if (tmp.compare("vtk_out_freq") == 0) {
				vtk_out_freq = getToken<cl_uint>(value);
			} else if (tmp.compare("vtk_out_name_base") == 0) {
				vtk_out_name_base = getToken<std::string>(value);
			} else if (tmp.compare("cl_plat_id") == 0) {
				cl_plat_id = getToken<cl_uint>(value);
			} else if (tmp.compare("cl_dev_id") == 0) {
				cl_dev_id = getToken<cl_uint>(value);
			} else if (tmp.compare("mass") == 0) {
				mass = getToken<cl_real>(value);
			} else if (tmp.compare("nu") == 0) {
				nu = getToken<cl_real>(value);
			} else if (tmp.compare("rho0") == 0) {
				rho0 = getToken<cl_real>(value);
			} else if (tmp.compare("k") == 0) {
				k = getToken<cl_real>(value);
			} else {
				std::cout << "Unknwon parameter:" << value;
			}


		}
		infile.close();

	} else {
		std::cout << "Error loading input file" << std::endl;
	}

}

void param::printInfo()
{
	std::cout << "/**************************************/" << std::endl;

	std::cout << "part_input_file: " << part_input_file << std::endl;
	std::cout << "dt: " << dt << std::endl;
	std::cout << "time_end: " << time_end << std::endl;
	std::cout << "part_out_freq: " << part_out_freq << std::endl;
	std::cout << "part_out_name_base: " << part_out_name_base << std::endl;
	std::cout << "vtk_out_freq: " << vtk_out_freq << std::endl;
	std::cout << "vtk_out_name_base: " << vtk_out_name_base << std::endl;
	std::cout << "cl_plat_id: " << cl_plat_id << std::endl;
	std::cout << "cl_dev_id: " << cl_dev_id << std::endl;
	std::cout << "XCells: " << XCells << std::endl;
	std::cout << "YCells: " << YCells << std::endl;
	std::cout << "ZCells: " << ZCells << std::endl;
	std::cout << "mass: " << mass << std::endl;
	std::cout << "nu: " << nu << std::endl;
	std::cout << "rho0: " << rho0 << std::endl;
	std::cout << "k: " << k << std::endl;

	std::cout << "/**************************************/" << std::endl;
}
