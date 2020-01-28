#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
int main()
{

	float Nx = 0.84;
	float Ny = 1.19;

	//float Nz = 1.5;
	float t=0.1;
	float k0 = 30000.;
	float rho0 = 7400.;
	float cp0 = 700.;
	float alpha0 = k0/(rho0*cp0); 
	float k=1000;
	float x0=0.5;
	float y0=1;
	float A=1/(4*M_PI*t*alpha0);
	int noPart=0;	

	std::ofstream   vtkFile("lessFluid.dat",std::ifstream::out);

	if (vtkFile.good())
    	{

		for (float i = 0.0; i <= 0.84; i+=0.015 ) 
		{
        		for (float j = 0.0; j <= 0.075;j+=0.015) 
			{
            		vtkFile << "4" << " " << i << "  " << j << "  " << 0.05 << " " << 0.0 << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			}
        	}

		for (float i = 0.0; i <= 0.075; i+=0.015 ) 
		{
        		for (float j = 0.09; j <= 2.19;j+=0.015) 
			{
            		vtkFile << "4" << " " << i << "  " << j << "  " << 0.05 << " " << 0.0 << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			}
        	}

		for (float i = 0.0; i <= 0.84; i+=0.015 ) 
		{
        		for (float j = 2.205; j <= 2.28;j+=0.015) 
			{
            		vtkFile << "4" << " " << i << "  " << j << "  " << 0.05 << " " << 0.0 << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			}
        	}

     	for (float i = 0.09; i <= Nx; i+=0.015 ) 
		{
        		for (float j = 0.09; j <= Ny;j+=0.015) 
			{
            		vtkFile << "1" << " " << i << "  " << j << "  " << 0.05 << " " << 0.0 << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			}
        	}


	for(float k = 0.855; k<=0.93; k+=0.015)
	{
			for(float l = 0.00; l<=2.28; l+=0.015)
			{
				
				vtkFile << "2" << " " << k << "  " << l << "  " << 0.05 << " " << "0" << " " << "0" << " " <<"0" << " " << 1500 << std::endl;
				noPart++;
			}
	}

	for(float k = 0.; k<=0.93; k+=0.015)
	{
			for(float l = 3; l<=3.51; l+=0.015)
			{
				
				vtkFile << "5" << " " << k << "  " << l << "  " << 0.05 << " " << "0" << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			}
	}

	/*for(float k = 0; k<2000; ++k)
	{
		
		vtkFile << "5" << " " << 0.5 << "  " << 3 << "  " << 0.05 << " " << "0" << " " << "0" << " " <<"0" << " " << 1500 << std::endl;
				noPart++;
	}*/
}

	/*for(float m = 3.025; m <= 5.0; m+=0.005)
	{
			//for(float n = 0.2; n<= 0.325; n+=0.025)
			//{	
				vtkFile << "2" << " " << 0.2 << "  " << m << "  " << 0.05 << " " << "0" << " " << "0" << " " <<"0" << " " << 2000 << std::endl;
				noPart++;
			//}
	}*/
	

	
	

	//vtkFile.seekp( std::ios_base::beg);
	std::cout << noPart ;
	//vtkFile << noPart << std::endl;
	vtkFile.close();

	return 0;
}
