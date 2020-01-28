#ifndef __PARAM_H__
#define	__RARAM_H__

#include <string>
#include <cstdlib>
#include <CL/cl.hpp>

#include "real.h"

// class to store simulation configuration
class param {
public:
	param(){};

	// possible by construct
//	explicit param( std::string pif, float dt, float tend, float ep, float sig, float lencell, cl_uint numX, cl_uint numY, cl_uint numZ, cl_uint pof, std::string ponb, 
//		cl_uint vof, std::string vonb, cl_uint pID, cl_uint devID, cl_uint groupsize);

	// load config file
	void	load_param(const std::string	filename);

	// print to console
	void	printInfo();

	~param(){}

private:
	param(const param& );
	void operator<<(const param& );

public:
	std::string 	part_input_file;
	cl_real		dt;
	cl_real		time_end;
	cl_uint		part_out_freq;
	std::string	part_out_name_base;
	cl_uint		vtk_out_freq;
	std::string	vtk_out_name_base;
	cl_uint		cl_plat_id;
	cl_uint		cl_dev_id;
	cl_uint		XCells;
	cl_uint		YCells;
	cl_uint		ZCells;
	cl_real		h;
	cl_real		mass;
	cl_real		nu;
	cl_real		rho0;
	cl_real		k;
};
#endif
