#ifndef __PARTICLE_H__
#define	__PARTICLE_H__

#include <string>
#include <vector>
#include <cstdlib>


#include <CL/cl.hpp>
#include "real.h"

// particle class
class particle {
public:
	particle();

	particle(cl_uint n);

	void	load_particles(const std::string	filename);

	~particle()
	{
		delete[] den;
		delete[] pres;
		delete[] tem;
		delete[] H;
		delete[] c_p;
		delete[] k;
		delete[] eps;
		delete[] pos;
		delete[] vel;
		delete[] force;
		delete[] pforce;
		delete[] vforce;
		delete[] id;
		
	}

	void erase();
	void resize(cl_uint n);

	// write format of md
	void writeMD(std::string filename);

	// VTK file write
	void writeVTK(std::string filename);

private:
	//forbid clonning
	particle(const particle& );
	void operator<<(const particle& );

public:
	cl_uint	N; // particle count

	cl_real*  den;
	cl_real*  pres;
	cl_real*  tem;
	cl_real* 	H;
	cl_real* 	c_p;
	cl_real* 	k;
	cl_real* 	eps;
	cl_real4* pos;
	cl_real4* vel; 
	cl_real4* force;
	cl_real4* pforce;
	cl_real4* vforce; 
	cl_uint*  id;

};
#endif
