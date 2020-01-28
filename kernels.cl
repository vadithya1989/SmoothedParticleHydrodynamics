
#include "real.h"

#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

__kernel void eulerIntegration
(
	__global const uint* id,
	__global real4* vel,
	__global real4* pos,
	__global const real4* force,
	const real mass,
	const real dt,
	const uint N,
	__global real* tem
)
{
	const uint i = get_global_id(0);
	if (i >= N) return;

	if (id[i] == 1 || id[i] == 3) {
		
		vel[i] += (force[i] / mass) * dt;	
		pos[i] += vel[i] * dt;

		//xtem[i] += (tem[i]) * dt;
		//printf("position --> %d %f %f %f %f %f %f\n", i, pos[i].x, pos[i].y, pos[i].z, vel[i].x, vel[i].y, vel[i].z);
	}
}
	
__kernel void blankCells
(
	__global int* Cells,
	const uint NCells
)
{
	const uint i = get_global_id(0);

	if (i >= NCells) return;

	Cells[i] = -1;
}

__kernel void sortParticles
(
	__global int* Cells,
	__global uint* Particles,
	__global const real4* pos,
	const uint XCells,
	const uint YCells,
	const uint ZCells,
	const real h,
	const uint N
)
{
	const uint i = get_global_id(0);
    if (i >= N) return;

	uint coord = pos[i].z / h;
    uint cID = coord * YCells;

    coord = pos[i].y / h;
    cID += coord;
    cID *= XCells;

    coord = pos[i].x / h;
    cID += coord;

    Particles[i] = atom_xchg(Cells + cID, i);

	//printf("particle %d is in %d\n", i, cID);
}

__kernel void calcForce
(
	__global const int* Cells,
	__global const int* Particles,
	__global const uint* cellsInDIM,
	__global const real4* pos,
	__global const real4* vel, 
	__global real4* force, 
	__global real4* pforce, 
	__global real4* vforce, 
	__global real* pres, 
	__global const real* den,
	const real mass, 
	const real h,
	//const real nu,
	__global real* tem,
	const real k,
//	const real rho0
	__global const uint* id
)
{
	int	cIndex = get_global_id(0);

	const uint XCells = cellsInDIM[0];
	const uint YCells = cellsInDIM[1];
	const uint ZCells = cellsInDIM[2];
	const uint NCells = XCells * YCells * ZCells;

	const real rho0 = 1900.f;
	
	if (cIndex >= NCells) return;

	int i = Cells[cIndex];	


	int cIndexX = cIndex % XCells;
	cIndex /= XCells;

	int cIndexY = cIndex % YCells;
	cIndex /= YCells;
	
	int cIndexZ = cIndex;

	int cIX = cIndexX;
	int cIY = cIndexY;
	int cIZ = cIndexZ;
	

	while (i != -1 ) {

		// "first compute pressure"
		//pres[i] = k * (den[i] - rho0);
		//real char_vel = 0.0011f;
		real ratio = den[i] / rho0;		
		pres[i] = k * rho0 / 7.f * ( pow(ratio, 7.f) - 1.f);

		// "re-calculate forces"
		pforce[i] = (real4)(0.f);
		vforce[i] = (real4)(0.f);
		force[i]  = (real4)(0.f);
		
	/*	real nu = 0.0f;

		if(tem[i] >= 1700.f && tem[i] < 2500.f )
		{
			nu = -0.000096667f * tem[i] + 0.194333f ;
		}
		else if(tem[i] <= 1700.f)
		{
			nu = 0.001f;
		}*/
		
		real nu = 0.1f;		
		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {

					cIndexX = cIX + dx;
					cIndexY = cIY + dy;
					cIndexZ = cIZ + dz;

					int cIndexCur = (cIndexZ * XCells * YCells) + (cIndexY * XCells) + cIndexX;
					if (cIndexCur >= 0 && cIndexCur < (XCells*YCells*ZCells ) ) {

						int j = Cells[cIndexCur];

						while (j != -1) {

							if (i != j)
							{
								real4 r = pos[i] - pos[j];
	
								real len_r = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

								if (len_r <= h && len_r > 0.f)
								{
									real h_five = h * h * h * h * h ;
									real M_PI = 3.141592654f;

									real pKernel = ( -30.f / ( M_PI * h_five ) ) * ((h - len_r) * (h - len_r)) * (1.f / len_r) * (M_PI * 4.f / 3.f * h * h * h);
									pforce[i] -= ((pres[i] + pres[j]) / (2.f * den[j])) * mass * pKernel * r;
									

									real vKernel = ( 20.0f / ( M_PI * h_five) ) * ( h - len_r );
									vforce[i] += nu * (vel[j] - vel[i]) * (mass / den[j]) * vKernel;

								}
							}
							j = Particles[j];
						}
					}
				}
			}
		}

		if (id[i] == 1 || id[i] == 3) {
			// "integreate all forces"
			real4 gHat = (real4)(0,-1,0,0); 					//modified
			real4 gforce = mass * gHat * 9.82f;
			real4 fringeForce = 0.0f;
			if(pos[i].y < 1.f && pos[i].y > 0.5f)
			{
			real v_0 = -0.35f;
			 fringeForce = gforce /(0.6f) * ((vel[i].y/v_0) - 1.0f);
			}
			force[i] = pforce[i] + vforce[i] + gforce ;//+ fringeForce;
			//printf("force --> %d %f %f %f %f %f %f %f %f %f %f %f\n", i, den[i], pres[i], pforce[i].x, pforce[i].y, pforce[i].z, vforce[i].x, vforce[i].y, vforce[i].z, pos[i].x, pos[i].y, pos[i].z);
		} else {
			force[i] = 0.f;
			pforce[i] = 0.f;
			vforce[i] = 0.f;
		}
		i = Particles[i];
	}
}

__kernel void calcTemp
(
	__global const int* Cells,
	__global const int* Particles,
	__global const uint* cellsInDIM,
	__global const real4* pos,
	__global const real* den,
	const real mass, 
	const real h,
	__global const uint* id,
	__global real* tem,
	__global real* H,
	__global real* c_p,
	__global real* k,
	__global real* eps
	
)

{

	
	int	cIndex = get_global_id(0);
	
	
	const uint XCells = cellsInDIM[0];
	const uint YCells = cellsInDIM[1];
	const uint ZCells = cellsInDIM[2];
	const uint NCells = XCells * YCells * ZCells;

	
	if (cIndex >= NCells) return;

	int i = Cells[cIndex];	
	
	if(id[i] == 1 || id[i] == 3)
	{
		int cIndexX = cIndex % XCells;
		cIndex /= XCells;

		int cIndexY = cIndex % YCells;
		cIndex /= YCells;
	
		int cIndexZ = cIndex;

		int cIX = cIndexX;
		int cIY = cIndexY;
		int cIZ = cIndexZ;
	

		while (i != -1 ) {
		
			//Calculate k
			real f1 = 0.0f, f2 = 0.0f;
			f1 = 43.079432359547056f - (0.003637526884958f * tem[i]) - (0.000052967611051f * (tem[i]*tem[i])) + (0.000000052645046f * (tem[i]*tem[i]*tem[i])) - (0.000000000013234f *(tem[i]*tem[i]*tem[i]*tem[i]));
			f2 = 4.114457824277824f +  (0.015038930375062f * tem[i]);
			k[i] = (tem[i]<=1800)*f1 + (tem[i]>1800)*f2;
		
			//Calculate c_p
			real g1 = 0.0f, g2 = 0.0f;
			g1 = 288.1686573077918f + (0.7408974274888f * tem[i]) - (0.0007665097216f * tem[i] * tem[i]) + (0.0000004230713f * tem[i] * tem[i] * tem[i]) - (0.0000000000854f * tem[i] * tem[i] * tem[i] * tem[i]);
			g2 =  818.2811987f;
			c_p[i] = (tem[i]<=1800)*g1 + (tem[i]>1800)*g2;

			//Calculate eps
			real h1 = 0.0f, h2 = 0.0f, h3 = 0.0f, h4 = 0.0f;
			h1 = 0.041365688318364f + (0.000112866817270f * tem[i]);
			h2 = -0.051666665902125f + (0.000208333332803f * tem[i]);
			h3 =  0.008469385594418f + (0.000136054423111f * tem[i]);
			h4 =  0.275f;
			eps[i] = (tem[i]<=1184)*h1 + (tem[i]>1184 && tem[i]<=1665)*h2 + (tem[i]>1665 && tem[i]<=1811)*h3 + (tem[i] > 1811)*h4;



			real tKernel = 0.0f;
			real t_sum_part = 0.0f;
			real hKernel = 0.0f;
			real h_sum_part = 0.0f;

			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					for (int dz = -1; dz <= 1; ++dz) {

						cIndexX = cIX + dx;
						cIndexY = cIY + dy;
						cIndexZ = cIZ + dz;

						int cIndexCur = (cIndexZ * XCells * YCells) + (cIndexY * XCells) + cIndexX;
						if (cIndexCur >= 0 && cIndexCur < (XCells*YCells*ZCells ) ) {

							int j = Cells[cIndexCur];

							while (j != -1) {

								if (i != j )
								{
									real4 r = pos[i] - pos[j];
	
									real len_r = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

									if (len_r <= h && len_r > 0.f)
									{
										real h_nine = h * h * h * h * h * h * h * h * h;
										real M_PI = 3.141592654f;
										real kk = ( 2.0968f / pos[i].y ) + 2.04255 ;
										real ccp = 7.0f;

										hKernel = -945.0f/(32.0f*M_PI*h_nine) * ((h * h) - (len_r * len_r)) * ((h * h) - (len_r * len_r));
										//h_sum_part += 0.001f *  mass / (den[i]*den[j]) * (4.0f * k[i] * k[j] / (k[i] + k[j] )) * (tem[i] - tem[j]) * hKernel ;
										h_sum_part += 0.001f * mass / (den[i]*den[j]) * (2.0f * kk) * (tem[i] - tem[j]) * hKernel ;

										//tKernel = -945.0f/(32.0f*M_PI*h_nine) * ((h * h) - (len_r * len_r)) * ((h * h) - (len_r * len_r));
										//t_sum_part += 0.001f * 1.0f/ccp * mass / (den[i]*den[j]) * (4.0f * kk * kk / (kk + kk )) * (tem[i] - tem[j]) * tKernel ;
									}
								}
								j = Particles[j];
							}
						}
					}
				}
			}
		
			//tem[i] += t_sum_part;
			H[i] += h_sum_part  ;
		
			float x1 = 1e2f * (3.672592527636285f + 0.000015688279392f * H[i]);
			float x2 =  1e2f * (3.485184453444313f + 0.000012220737794f * H[i]);
			float x3 = 1811.f;
			tem[i] =  (0.f <= H[i] && H[i] <= 944393.1f)*x1 + (944393.1f < H[i] && H[i] <= 1196721.6f)*x3 + ( 1196721.6f < H[i])*x2;
			i = Particles[i];
	 
		}
	}
}




__kernel void calcDensity
(
	__global const int* Cells,
	__global const int* Particles,
	__global const uint* cellsInDIM,
	__global const real4* pos, 
	__global real* den, 
	const real mass,
	const real h,
	__global const uint* id
)
{
	int	cIndex = get_global_id(0);

	const uint XCells = cellsInDIM[0];
	const uint YCells = cellsInDIM[1];
	const uint ZCells = cellsInDIM[2];
	const uint NCells = XCells * YCells * ZCells;

	if (cIndex >= NCells) return;

	int i = Cells[cIndex];

	int cIndexX = cIndex % XCells;
	cIndex /= XCells;

	int cIndexY = cIndex % YCells;
	cIndex /= YCells;
	
	int cIndexZ = cIndex;

	int cIX = cIndexX;
	int cIY = cIndexY;
	int cIZ = cIndexZ;
	
	

	while (i != -1) {

		// "re-calculate density"
		den[i] = (real)0.f;
		real dKernel = 0.f;

		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {

					cIndexX = cIX + dx;
					cIndexY = cIY + dy;
					cIndexZ = cIZ + dz;

					int cIndexCur = (cIndexZ * XCells * YCells) + (cIndexY * XCells) + cIndexX;
		//		printf("id-> %d\n",cIndexCur);
					if (cIndexCur >= 0 && cIndexCur < (XCells*YCells*ZCells ) ) {

						int j = Cells[cIndexCur];
						while (j != -1) {

							real4 r = pos[i] - pos[j];
//				printf("pos--> %d %f %f %f %f %f %f %f %f %f\n", i, r.x, r.y, r.z, pos[i].x, pos[i].y, pos[i].z, pos[j].x, pos[j].y, pos[j].z);
							real len_r2 = (r.x*r.x + r.y*r.y + r.z*r.z);
							real M_PI = 3.141592654f;

							if (len_r2 <= (h*h) && len_r2 > 0.f)
							{
								real len_r = sqrt(len_r2);	
								real h_eight =  h * h * h * h * h * h * h * h;
								dKernel += (4.f / ( M_PI * h_eight)) * ((h * h) - (len_r * len_r)) * ((h * h) - (len_r * len_r)) * ((h * h) - (len_r 							* len_r));
							}
							j = Particles[j];
						}
					}
				}
				
			}
		}
		//if (id[i] == 1) {
			den[i] += mass * dKernel;
		//} else {
			//den[i] = 1752.f;
			//den[i] = 205.f;
		//}
		i = Particles[i];
	}
}

__kernel void applyPeriodic
(
    __global const real4* pos,
    __global real4* vel,
    unsigned int N,
	__global real* tem,
	__global uint* id,
	__global real4* vforce,
	__global real4* pforce,
	__global real* H,
	__global uint* dynID
)
{
    unsigned int i = get_global_id(0);

     if (i >= N) return;
	//printf("count ---> %d\n", count);
	
    /*if ( pos[i].x < 0.) {
        vel[i].x = fabs(vel[i].x);
    }

    if ( pos[i].x >= 3 ) {
        vel[i].x = -fabs(vel[i].x);
    }

    if ( pos[i].y < 0.) {
        vel[i].y = fabs(vel[i].y);
    }

    if ( pos[i].y >= 3 ) {
        vel[i].y = -fabs(vel[i].y);
    }

    if ( pos[i].z < 0.) {
        vel[i].z = fabs(vel[i].z);
    }

    if ( pos[i].z >= 1.025 ) {
        vel[i].z = -fabs(vel[i].z);
    }*/
	

	
	

	if(pos[i].y >= 0.075f && pos[i].y <= 0.2f)
	{
		vel[i].y = -0.07f;
	}		
	

	if ( pos[i].y < 0.075f && id[i] == 1) {
		
		id[i] = 5;
		tem[i] = 2000.f;
		real f1 = 0.0f, f2 = 0.0f;
		f1 = ( -2.274393289788107f +  0.006303449978286f * tem[i]) * 1e5f;
		f2 = (-2.851860920398242f + 0.008182812010365f * tem[i] ) * 1e5f; 
		H[i] = (tem[i]<=1800.f && tem[i] >= 298.f)*f1 + (tem[i]>1800.f)*f2;
		
		pos[i].y = 5.f;
		pos[i].x = 0.5f ;
		/*vel[i].x = 0.0f;
		vel[i].y = 0.0f;	
		vforce[i].x = 0.0f;vforce[i].y = 0.0f;vforce[i].z = 0.0f;
		pforce[i].x = 0.0f;pforce[i].y = 0.0f;pforce[i].z = 0.0f;*/
		
    }
	
}


__kernel void applySpring
(
	__global real4* pos, 
	__global real4* force, 
	__global const real4* vel, 
	const unsigned int XCells, 
	const unsigned int YCells, 
	const unsigned int ZCells, 
	const real  h,
	unsigned int N,
	__global const uint* id
)
{
    unsigned int i = get_global_id(0);

	
    if (i >= N) return;

        real4 dij;

		real xMax = 0.855f; 				//modified
		real yMax = 5.0f;				//modified
		real zMax = 0.4f;

        dij.x = min( (xMax - pos[i].x ), pos[i].x - 0.075f);
        dij.y = min( (yMax - pos[i].y ), pos[i].y - 0.075f);
	//dij.y = (yMax-pos[i].y);
        dij.z = min( (zMax - pos[i].z ), pos[i].z);

        real4 repul; 
        repul.x = 0.f; 
        repul.y = 0.f; 
        repul.z = 0.f;

		// "define constants"
		real k_s = 300.f;
		real k_d = 50.f;

        	if (dij.x < (0.015f) ) {
            repul.x     = 2.f * k_s * ( h - dij.x) + k_d * -2.f * ( h * vel[i].x ) ;
            if (dij.x == pos[i].x) {
                force[i].x += fabs(repul.x);
            } else {
                force[i].x -= fabs(repul.x);
            }
        	}
      	
		if (dij.y < (0.015f)) {
            repul.y     = 2.f * k_s * ( h - dij.y) + k_d * -2.f * ( h * vel[i].y ) ;
            if (dij.y == pos[i].y) {
                force[i].y += fabs(repul.y);
            } else {
                force[i].y -= fabs(repul.y);
            }
        	}
	

        if (dij.z < (0.015f)) {
            repul.z     = 2.f * k_s * ( h - dij.z) + k_d * -2.f * ( h * vel[i].z ) ;
            if (dij.z == pos[i].z) {
                force[i].z += fabs(repul.z);
            } else {
                force[i].z -= fabs(repul.z);
            }
        }

}


__kernel void calcSurfaceForce
(
	__global const int* Cells,
	__global const int* Particles,
	__global const uint* cellsInDIM,
	__global const real4* pos,
	__global const real4* vel, 
	__global real4* force, 
	__global real4* pforce, 
	__global real4* vforce, 
	__global real* pres, 
	__global const real* den,
	const real mass, 
	const real h,
	const real nu,
	const real k,
//	const real rho0
	__global const uint* id
)
{
	int	cIndex = get_global_id(0);

	const uint XCells = cellsInDIM[0];
	const uint YCells = cellsInDIM[1];
	const uint ZCells = cellsInDIM[2];
	const uint NCells = XCells * YCells * ZCells;

	if (cIndex >= NCells) return;

	int i = Cells[cIndex];	

	int cIndexX = cIndex % XCells;
	cIndex /= XCells;

	int cIndexY = cIndex % YCells;
	cIndex /= YCells;
	
	int cIndexZ = cIndex;

	int cIX = cIndexX;
	int cIY = cIndexY;
	int cIZ = cIndexZ;
	

	while (i != -1 ) {

		// "re-calculate forces"
		real4 ni = (real4)0.f;
		real4 lapci = (real4)0.f;

		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {

					cIndexX = cIX + dx;
					cIndexY = cIY + dy;
					cIndexZ = cIZ + dz;

					int cIndexCur = (cIndexZ * XCells * YCells) + (cIndexY * XCells) + cIndexX;
					if (cIndexCur >= 0 && cIndexCur < (XCells*YCells*ZCells ) ) {

						int j = Cells[cIndexCur];

						while (j != -1) {

							real4 r = pos[i] - pos[j];

							real len_r = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

							if (len_r <= h)
							{
								real h_nine = h * h * h * h * h * h * h * h * h;
								real M_PI = 3.141592654f;

								real sKernel = ( -945.f / ( 32.f * M_PI * h_nine ) ) *  ((h * h) - (len_r * len_r)) * ((h * h) - (len_r * len_r));
								real temp = ( -945.f / ( 32.f * M_PI * h_nine ) ) *  ((h * h) - (len_r * len_r)) * ((3.f* h*h) - (7.f* len_r * len_r));

								ni += (mass / den[j]) * sKernel * r;
								lapci += (mass / den[j]) * temp;
							}
							j = Particles[j];
						}
					}
				}
			}
		}

		real len_ni = sqrt( ni.x*ni.x + ni.y*ni.y + ni.z*ni.z );

		// "paramter l"
		real l = 7.f;

		if ( len_ni > l) {
			// "integreate surface tension force"
			real sigma = 0.0000728f;
			real4 sforce = (-sigma * lapci * (ni / len_ni));
//		printf("sforface--> %d %f %f %f\n",i, sforce.x, sforce.y, sforce.z);
			force[i] += sforce;
			
			
		}
		i = Particles[i];
	}
}


__kernel void dynReservoir
(
	__global uint* id,
	__global uint* count,
	__global uint* dynID
)
{

	unsigned int i = get_global_id(0);

	if(id[i] == 5)
	{
		count[0] = count[0] + 1;
		dynID[count[0] - count[0]] = i;
		//printf("%d\n", i);
	}
	
}


__kernel void dynReservoirOne
(
    __global const real4* pos,
    __global real4* vel,
    unsigned int N,
	__global real* tem,
	__global uint* id,
	__global real4* vforce,
	__global real4* pforce,
	__global real* H,
	__global uint* dynID
)
{
    unsigned int i = get_global_id(0);

     if (i >= N) return;


	if(i == dynID[0] && pos[i].y < 4.5f)
	{
		id[i] = 1;
		pos[i].x = 0.15f;
		pos[i].y = 0.75f;
		vel[i].y = -0.0f;
		vel[i].x = 1.0f;
		vforce[i].x = 0.0f;vforce[i].y = 0.0f;vforce[i].z = 0.0f;
		pforce[i].x = 0.0f;pforce[i].y = 0.0f;pforce[i].z = 0.0f;
		//printf("adfg%d\n", i);
		
	}

}


__kernel void movingWall
(
    __global real4* pos,
	__global uint* id
)
{
	
	unsigned int i = get_global_id(0);

	if(id[i] == 2)
	{
				
		pos[i].y = pos[i].y + 1.015;
		id[i] = 6;
		//printf("df\n");
	}

	/*if(id[i] == 6)
	{
		id[i] = 2;		
		pos[i].y = pos[i].y - 1.015;
		printf("df\n");
		
	}*/

}


