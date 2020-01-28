#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <limits>
#include <boost/lexical_cast.hpp>

template<typename T>
T getToken(const std::string& line)
{
	T	value;
    std::stringstream    ss;
	ss << line;
    ss >> value;

	return value;
}

particle::particle()
:
N(0),
den(NULL),
pres(NULL),
tem(NULL),
H(NULL),
c_p(NULL),
k(NULL),
eps(NULL),
pos(NULL),
vel(NULL),
force(NULL),
pforce(NULL),
vforce(NULL),
id(0)
{}

particle::particle(cl_uint n)
:
N(n)
{

	resize(n);
}

void particle::resize(cl_uint n)
{
	N = n;

	den 		= new cl_real[N];
	pres 	= new cl_real[N];
	tem 		= new cl_real[N];
	H   		= new cl_real[N];
	c_p  	= new cl_real[N];
	k   		= new cl_real[N];
	eps  	= new cl_real[N];
	pos 		= new cl_real4[N];
	vel 		= new cl_real4[N];
	force 	= new cl_real4[N];
	pforce 	= new cl_real4[N];
	vforce 	= new cl_real4[N];
	id 		= new cl_uint[N];

	for (unsigned int i = 0; i < N; ++i)
	{
		den[i] = 0;
		pres[i] = 0;
		tem[i] = 0;
		H[i] = 0;
		c_p[i] = 0;k[i] = 0;eps[i] = 0;
		pos[i].x = pos[i].y = pos[i].z = 0;
		vel[i].x = vel[i].y = vel[i].z = 0;
		force[i].x = force[i].y = force[i].z = 0;
		pforce[i].x = pforce[i].y = pforce[i].z = 0;
		vforce[i].x = vforce[i].y = vforce[i].z = 0;
		id[i] = 0;
	}

}

void particle::erase()
{
	delete [] den;
	delete [] pres;
	delete [] tem;
	delete [] H;
	delete [] c_p;
	delete [] k;
	delete [] eps;
	delete [] pos;
	delete [] vel;
	delete [] force;
	delete [] pforce;
	delete [] vforce;
	delete [] id;

	id = 0; 
	den = 0; pres = 0; tem = 0; H = 0; pos = 0; vel = 0; force = 0; pforce=0; vforce=0;c_p = 0;k = 0;eps = 0;
}

void    particle::load_particles(const std::string  filename)
{
    std::ifstream   infile(filename.c_str(),std::ifstream::in);
    if (infile.good())
    {
        std::string line, tmp;
        getline(infile, line);
	    std::stringstream	s1(line);
	    s1 >> N;
        if (infile) {
            erase();
			resize(N);
            for (unsigned int i = 0; i < N; ++i)
            {
                	getline(infile, line);
        			std::stringstream    s2(line);
				den[i] = 0;
				pres[i] = 0;
				H[i] = 0;
				c_p[i] = 0;
				k[i] = 0;
				eps[i] = 0;
		        	s2 >> id[i] >> pos[i].x >> pos[i].y >> pos[i].z >> vel[i].x >> vel[i].y >> vel[i].z >> tem[i] ; //modified
				
				double f1 = 0.0, f2 = 0.0;
				f1 = ( -2.274393289788107 +  0.006303449978286 * tem[i]) * 1e5;
				f2 = (-2.851860920398242 + 0.008182812010365 * tem[i] ) * 1e5; 
				H[i] = (tem[i]<=1800 && tem[i] >= 298)*f1 + (tem[i]>1800)*f2; 
				force[i].x = force[i].y = force[i].z = 0;
				pforce[i].x = pforce[i].y = pforce[i].z = 0;
				vforce[i].x = vforce[i].y = vforce[i].z = 0;
            }
        } else {
			std::cerr << "Error in particles file" << std::endl;
		}

    } else {
        std::cout << "Error loading input file" << std::endl;
			exit(-1);
    }
}

void particle::writeMD(std::string filename)
{
    std::ofstream   mdfile(filename.c_str(),std::ifstream::out);
    if (mdfile.good())
    {
		mdfile << N << std::endl;

		for (unsigned int i = 0; i < N; ++i)
		{
			mdfile << " " << den[i] << " " << pres[i] << " " << tem[i] << " " << pos[i].x << " " << pos[i].y << " " << pos[i].z << " " 
					<< vel[i].x << " " << vel[i].y << " " << vel[i].z << " "
					<< force[i].x << " " << force[i].y << " " << force[i].z << " " 
					<< pforce[i].x << " " << pforce[i].y << " " << pforce[i].z << " "
					<< vforce[i].x << " " << vforce[i].y << " " << vforce[i].z << std::endl;

		}
	}
	mdfile.close();
}

void particle::writeVTK(std::string filename)
{
    std::ofstream   vtkFile(filename.c_str(),std::ifstream::out);
    if (vtkFile.good())
    {
		vtkFile.precision(15);

		// Header
		vtkFile << "# vtk DataFile Version 3.0" << std::endl;
		vtkFile << "SPH particles" << std::endl;
		vtkFile << "ASCII" << std::endl;
		vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		// write grid from position
		vtkFile << "POINTS " << size_t(N) << " float" << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << std::fixed << "	" << pos[i].x << "	" << pos[i].y << "	" << pos[i].z << std::endl;
		}

		// write density
		vtkFile << "POINT_DATA " << size_t(N) << std::endl;
		vtkFile << "SCALARS density float 1"  << std::endl;
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << den[i] << std::endl;
		}

		//temperature
		vtkFile << "SCALARS temperature float 1"  << std::endl;   	//modified
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << tem[i] << std::endl;
		}

		//Enthalpy
		vtkFile << "SCALARS enthalpy float 1"  << std::endl;   	//modified
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << H[i] << std::endl;
		}

		// write ID
		vtkFile << "SCALARS id float 1"  << std::endl;
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << id[i] << std::endl;
		}

		// write pressure
		vtkFile << "SCALARS pressure float 1"  << std::endl;
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << pres[i] << std::endl;
		}

		/*// write temperatur
		vtkFile << "SCALARS temperatur float 1"  << std::endl;
		vtkFile << "LOOKUP_TABLE default"  << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << tem[i] << std::endl;
		}*/

		// write velocity
		vtkFile << "VECTORS velocity float" << std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << vel[i].x << "	" << vel[i].y << "	" << vel[i].z << std::endl;
		}

		// write force
		vtkFile << "VECTORS force float"<< std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << force[i].x << "	" << force[i].y << "	" << force[i].z << std::endl;
		}

		// write pressure force
		vtkFile << "VECTORS pforce float"<< std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << pforce[i].x << "	" << pforce[i].y << "	" << pforce[i].z << std::endl;
		}

		// write viscous force
		vtkFile << "VECTORS vforce float"<< std::endl;
		for (unsigned int i = 0; i < N; ++i )
		{
			vtkFile << "	" << vforce[i].x << "	" << vforce[i].y << "	" << vforce[i].z << std::endl;
		}
	}

	vtkFile.close();
}
