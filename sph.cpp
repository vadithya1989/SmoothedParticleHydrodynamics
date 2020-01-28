
#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <fstream>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include "particle.h"
#include "param.h"
#include "real.h"

#include "useOpenCL.h"

// help function for writing results 
void writeOutput( size_t step, const param& prm, particle& md, std::string method )
{
	if ( method.compare("MD") == 0)
	{
        std::stringstream s;
        std::string file_MD;
        s << prm.part_out_name_base << step << ".md";
        s >> file_MD;
        std::cout << "...saving_MD";
        md.writeMD(file_MD);
	} else if ( method.compare("VTK") == 0) {
        std::stringstream s;
        std::string file_VTK;
        s << prm.vtk_out_name_base << step << ".vtk";
        s >> file_VTK;
        std::cout << "...saving_VTK";
        md.writeVTK(file_VTK);
	} else {
		std::cout << "Error: Unknown write format" << std::endl;
		exit(-1);
	}
}


int main(int argc, char* argv[])
{

	try
	{

	if (argc < 2) {
		std::cout << "./sph <cfg_file>" << std::endl;
		return -1;
	}

	// init config
	param prm;
	prm.load_param(argv[1]);
	prm.printInfo();

	// init particle
	particle md;
	md.load_particles(prm.part_input_file);

	// handy OpenCL class
	OpenCL::OpenCL* ocl = &OpenCL::OpenCL::instance();
	ocl->initById(prm.cl_plat_id, true, prm.cl_dev_id);
	std::vector<cl::Device> devices = ocl->devices();
	cl::Context context = ocl->context();
	cl::CommandQueue cmdQ = ocl->queue();
	
	// cell Linked List
	cl_uint NCells = prm.XCells * prm.YCells * prm.ZCells;
	cl_int Cells[NCells];

	// particle Linked List
	cl_uint Particles[md.N];

	// initalize cells
	for (cl_uint i = 0; i < NCells; ++i)
	{
		Cells[i] = -1;
	}

	// work around for parameters
	cl_uint* cellsInDIM = new cl_uint[3];
	cellsInDIM[0] = prm.XCells;
	cellsInDIM[1] = prm.YCells;
	cellsInDIM[2] = prm.ZCells;

	//**********************************FOR DYNAMIC RESERVOIR************************************//

	cl_uint* count = new cl_uint[1];
	count[0] = 0;

	cl_uint* dynID = new cl_uint[6000];
	for(cl_uint i = 0; i < 6000 ; ++i)
	{
		dynID[i] = 0;
	}

	cl_uint* dynID1 = new cl_uint[6000];
	for(cl_uint i = 0; i < 6000 ; ++i)
	{
		dynID1[i] = 0;
	}

	// buffer for count and dynID
	cl::Buffer	buf_count(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * 1);
	cl::Buffer	buf_dynID(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * 6000);

	// write to the above buffers
	cmdQ.enqueueWriteBuffer(buf_count, CL_TRUE, 0, sizeof(cl_uint) * 1, &count[0]);
	cmdQ.enqueueWriteBuffer(buf_dynID, CL_TRUE, 0, sizeof(cl_uint) * 6000, &dynID[0]);

	//*******************************************************************************************//

	// create cl::Buffer
	cl::Buffer	buf_den(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_pres(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_tem(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_H(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_c_p(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_k(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_eps(context, CL_MEM_READ_WRITE, sizeof(cl_real) * md.N);
	cl::Buffer	buf_pos(context, CL_MEM_READ_WRITE, sizeof(cl_real4) * md.N);
	cl::Buffer	buf_vel(context, CL_MEM_READ_WRITE, sizeof(cl_real4) * md.N);
	cl::Buffer	buf_force(context, CL_MEM_READ_WRITE, sizeof(cl_real4) * md.N);
	cl::Buffer	buf_pforce(context, CL_MEM_READ_WRITE, sizeof(cl_real4) * md.N);
	cl::Buffer	buf_vforce(context, CL_MEM_READ_WRITE, sizeof(cl_real4) * md.N);
	cl::Buffer	buf_id(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * md.N);

	cl::Buffer	buf_Cells(context, CL_MEM_READ_WRITE, sizeof(cl_int) * NCells);
	cl::Buffer	buf_Particles(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * md.N);
	cl::Buffer	buf_cellsInDIM(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * 3);

	// write buffer
	cmdQ.enqueueWriteBuffer(buf_den, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.den[0]);
	cmdQ.enqueueWriteBuffer(buf_pres, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.pres[0]);
	cmdQ.enqueueWriteBuffer(buf_tem, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.tem[0]);
	cmdQ.enqueueWriteBuffer(buf_H, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.H[0]);
	cmdQ.enqueueWriteBuffer(buf_c_p, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.c_p[0]);
	cmdQ.enqueueWriteBuffer(buf_k, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.k[0]);
	cmdQ.enqueueWriteBuffer(buf_eps, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.eps[0]);
	cmdQ.enqueueWriteBuffer(buf_pos, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.pos[0]);
	cmdQ.enqueueWriteBuffer(buf_vel, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.vel[0]);
	cmdQ.enqueueWriteBuffer(buf_force, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.force[0]);
	cmdQ.enqueueWriteBuffer(buf_pforce, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.pforce[0]);
	cmdQ.enqueueWriteBuffer(buf_vforce, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.vforce[0]);
	cmdQ.enqueueWriteBuffer(buf_id, CL_TRUE, 0, sizeof(cl_uint) * md.N, &md.id[0]);

	cmdQ.enqueueWriteBuffer(buf_Cells, CL_TRUE, 0, sizeof(cl_int) * NCells, &Cells[0]);
	cmdQ.enqueueWriteBuffer(buf_Particles, CL_TRUE, 0, sizeof(cl_uint) * md.N, &Particles[0]);
	cmdQ.enqueueWriteBuffer(buf_cellsInDIM, CL_TRUE, 0, sizeof(cl_int) * 3, &cellsInDIM[0]);

	// kernel
	std::ifstream kin("kernels.cl");

    std::istreambuf_iterator<char> begin(kin), end;
    std::string kernels(begin, end);

	cl::Program::Sources source;

	source.push_back(std::make_pair(kernels.c_str(), kernels.length()+1));

	cl::Program program(context, source);

	try {
    	program.build(devices,"-I. -g");
    } catch (cl::Error& err) {
        std::cerr << "Building failed, " << err.what() << "(" << err.err() << ")"
              << "\nRetrieving build log\n"
              << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[prm.cl_dev_id])
              << "\n";
	}

	cl::Kernel eulerIntegration(program,"eulerIntegration");
	cl::Kernel calcForce(program,"calcForce");
	cl::Kernel calcTemp(program,"calcTemp");
	cl::Kernel calcDensity(program,"calcDensity");
	cl::Kernel blankCells(program,"blankCells");
	cl::Kernel sortParticles(program,"sortParticles");
	cl::Kernel applyPeriodic(program,"applyPeriodic");
	cl::Kernel applySpring(program,"applySpring");
	cl::Kernel calcSurfaceForce(program,"calcSurfaceForce");
	cl::Kernel dynReservoir(program,"dynReservoir");
	cl::Kernel dynReservoirOne(program,"dynReservoirOne");
	cl::Kernel movingWall(program,"movingWall");

	cl::KernelFunctor funcEuler		( eulerIntegration, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcForce		( calcForce, 	cmdQ, cl::NullRange, cl::NDRange(NCells), cl::NDRange(1) );
	cl::KernelFunctor funcTemp		( calcTemp, 	cmdQ, cl::NullRange, cl::NDRange(NCells), cl::NDRange(1) );
	cl::KernelFunctor funcDensity		( calcDensity, cmdQ, cl::NullRange, cl::NDRange(NCells), cl::NDRange(1) );
	cl::KernelFunctor funcBlankCells	( blankCells,   cmdQ, cl::NullRange, cl::NDRange(NCells), cl::NDRange(1) );
	cl::KernelFunctor funcSortParticles( sortParticles, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcPeriodic	( applyPeriodic, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcSpring		( applySpring, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcSurface		( calcSurfaceForce, cmdQ, cl::NullRange, cl::NDRange(NCells), cl::NDRange(1) );
	cl::KernelFunctor funcReservoir	( dynReservoir, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcReservoirOne	( dynReservoirOne, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );
	cl::KernelFunctor funcMovingWall	( movingWall, cmdQ, cl::NullRange, cl::NDRange(md.N), cl::NDRange(1) );

	// performance count
	cl::Event event;

	//std::cout << "clear Cells..." << std::endl;
    	event = funcBlankCells(buf_Cells, NCells);
    	event.wait();

    	//std::cout << "sorting..." << std::endl;
    	event = funcSortParticles(buf_Cells, buf_Particles, buf_pos, prm.XCells, prm.YCells, prm.ZCells, prm.h, md.N);
    	event.wait();

	//std::cout << "density..." << std::endl;
	
	event = funcDensity(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos, buf_den, prm.mass, prm.h, buf_id); 
    	event.wait();

	//std::cout << "temp..." << std::endl;
	event = funcTemp(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos,  buf_den, prm.mass, prm.h, buf_id, buf_tem, buf_H, buf_c_p, buf_k, buf_eps); 
     event.wait();

	//std::cout << "force..." << std::endl;
	event = funcForce(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos, buf_vel, buf_force, buf_pforce, buf_vforce, buf_pres, buf_den, prm.mass, prm.h ,/*prm.nu,*/ buf_tem, prm.k, /*prm.rho0*/ buf_id); 
	event.wait();


	int  step = 0, mdCount = 0, vtkCount = 0;

	// "time loop"
	for (cl_real t = 0; t <= prm.time_end; t += prm.dt, ++step)
	{
		

		std::cout << "\nTime = " << t ;

		//std::cout << "euler..." << std::endl;
		event = funcEuler(buf_id, buf_vel, buf_pos, buf_force, prm.mass, prm.dt, md.N, buf_tem);
		event.wait();

		//std::cout << "clear Cells..." << std::endl;
    		event = funcBlankCells(buf_Cells, NCells);
    		event.wait();
	
    		//std::cout << "sorting..." << std::endl;
    		event = funcSortParticles(buf_Cells, buf_Particles, buf_pos, prm.XCells, prm.YCells, prm.ZCells, prm.h, md.N);
    		event.wait();

		//std::cout << "density..." << std::endl;
		event = funcDensity(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos, buf_den, prm.mass, prm.h, buf_id); 
    		event.wait();

		//std::cout << "force..." << std::endl;
		event = funcForce(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos, buf_vel, buf_force, buf_pforce, buf_vforce, buf_pres, buf_den, prm.mass, prm.h, /*prm.nu */buf_tem , prm.k, /*prm.rho0*/ buf_id);		
		event.wait();

		//std::cout << "temp..." << std::endl;
		event = funcTemp(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos,  buf_den, prm.mass, prm.h, buf_id, buf_tem, buf_H, buf_c_p, buf_k, buf_eps); 
     	event.wait();

		//std::cout << "force..." << std::endl;
		event = funcSurface(buf_Cells, buf_Particles, buf_cellsInDIM, buf_pos, buf_vel, buf_force, buf_pforce, buf_vforce, buf_pres, buf_den, prm.mass, prm.h, prm.nu, prm.k, /*prm.rho0*/ buf_id);
    		event.wait();


		if(step % 5 ==0)
		{
			//std::cout << "dyn..." << std::endl;
			event = funcReservoir(buf_id, buf_count, buf_dynID);
	    		event.wait();

			//cmdQ.enqueueReadBuffer(buf_dynID, CL_TRUE, 0, sizeof(cl_uint) * 6000, &dynID1[0]);
			//std::cout << std::endl << dynID1[0] << std::endl;

			//std::cout << "dyn1..." << std::endl;
			event = funcReservoirOne(buf_pos, buf_vel, md.N, buf_tem, buf_id, buf_vforce, buf_pforce, buf_H, buf_dynID);
		    	event.wait();
		
		}

		//std::cout << "periodic..." << std::endl;
		event = funcPeriodic(buf_pos, buf_vel, md.N, buf_tem, buf_id, buf_vforce, buf_pforce, buf_H, buf_dynID);
	    	event.wait();

		//std::cout << "spring bc..." << std::endl;
		event = funcSpring(buf_pos, buf_force, buf_vel, prm.XCells, prm.YCells, prm.ZCells, prm.h, md.N, buf_id);
		event.wait();

		/*if(step % 20 == 0)
		{
			//std::cout << "spring bc..." << std::endl;
			event = funcMovingWall(buf_pos, buf_id);
			event.wait();
		}*/

		

		/*std::vector< cl_double > enthalpy(md.N, 0.0);
		cmdQ.enqueueReadBuffer(buf_tem, CL_TRUE, 0, sizeof(cl_real) * md.N, &enthalpy[0]);
		//calcEnthalpy(enthalpy);
		std::cout  << std::endl << enthalpy.at(50) << std::endl;*/

		// read whenever
		if ((step % prm.part_out_freq == 0) ||  (step % prm.vtk_out_freq == 0))
		{
			cmdQ.enqueueReadBuffer(buf_den, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.den[0]);
			cmdQ.enqueueReadBuffer(buf_pres, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.pres[0]);
			cmdQ.enqueueReadBuffer(buf_tem, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.tem[0]);
			cmdQ.enqueueReadBuffer(buf_id, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.id[0]);
			cmdQ.enqueueReadBuffer(buf_H, CL_TRUE, 0, sizeof(cl_real) * md.N, &md.H[0]);
			cmdQ.enqueueReadBuffer(buf_pos, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.pos[0]);
			cmdQ.enqueueReadBuffer(buf_vel, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.vel[0]);
			cmdQ.enqueueReadBuffer(buf_force, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.force[0]);
			cmdQ.enqueueReadBuffer(buf_pforce, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.pforce[0]);
			cmdQ.enqueueReadBuffer(buf_vforce, CL_TRUE, 0, sizeof(cl_real4) * md.N, &md.vforce[0]);
		}


		// write md file
		if ( step >= 0 && (step % prm.part_out_freq == 0))
		{
			writeOutput(mdCount, prm, md, "MD");
			++mdCount;
		}

		// write vtk file
		if ( step >= 0 && (step % prm.vtk_out_freq == 0))
		{
			writeOutput(vtkCount, prm, md, "VTK");
			++vtkCount;
		}

	}

	std::cout << "\n\nFinished." << std::endl;
	} catch(cl::Error& err){
		 std::cerr << "An OpenCL error occured, " << err.what()
                << "\nError num of " << err.err() << "\n";
		return -1;
    }

	return 0;  
}
