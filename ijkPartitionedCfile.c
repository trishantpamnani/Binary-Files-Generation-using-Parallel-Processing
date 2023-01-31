#include "mpi.h"
using namespace std;
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include "TECIO.h"
#define pi acos(-1)

#ifndef NULL
#define NULL 0
#endif

#define XDIM 100
#define YDIM 100
#define ZDIM 100




INTEGER4 outputVarData(
    double* var,
    INTEGER4 iDim, INTEGER4 jDim, INTEGER4 kDim,
    INTEGER4 iMin, INTEGER4 jMin, INTEGER4 kMin,
    INTEGER4 iMax, INTEGER4 jMax, INTEGER4 kMax);

int main(int argc, char** argv)
{
    INTEGER4 Debug = 1;
    INTEGER4 VIsDouble = 1;
    INTEGER4 FileType = 0;
    INTEGER4 FileFormat = 1; // SZPLT supported for partitioned zones
    INTEGER4 I = 0; // Used to track return codes
    
    //MPI initiation
    MPI_Init(&argc, &argv);
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    int commSize;
    MPI_Comm_size(mpiComm, &commSize);
    int commRank;
    MPI_Comm_rank(mpiComm, &commRank);
    
    

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"IJK Ordered Zone",
        (char*)"X Y Z rho ux uy uz",
        (char*)"MPIijkpartitioned.szplt",
        (char*)".",
        &FileFormat,
        &FileType,
        &Debug,
        &VIsDouble);

    INTEGER4 mainRank = 0;
    I = TECMPIINIT142(&mpiComm, &mainRank);


    //declaration of input arrays
    double* x = new double[XDIM * YDIM * ZDIM];
    double* y = new double[XDIM * YDIM * ZDIM];
    double* z = new double[XDIM * YDIM * ZDIM];
    double* rho = new double[XDIM * YDIM * ZDIM];
    double* ux = new double[XDIM * YDIM * ZDIM];
    double* uy = new double[XDIM * YDIM * ZDIM];
    double* uz = new double[XDIM * YDIM * ZDIM];

    //parameters for Teczne142
    INTEGER4 IMax = XDIM;
    INTEGER4 JMax = YDIM;
    INTEGER4 KMax = ZDIM;
    INTEGER4 ICellMax = 0;
    INTEGER4 JCellMax = 0;
    INTEGER4 KCellMax = 0;
    INTEGER4 DIsDouble = 1;
    double   SolTime = 360.0;
    INTEGER4 StrandID = 0;      /* StaticZone */
    INTEGER4 unused = 0;      // ParentZone is no longer used
    INTEGER4 IsBlock = 1;      /* Block */
    INTEGER4 NFConns = 0;
    INTEGER4 FNMode = 0;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 valueLocations[] = { 1, 1, 1, 1,1,1,1 };
    INTEGER4 ShrConn = 0;
    double U0 = 0.05;

    for (int k = 0; k < ZDIM; k++)
    {

        for (int j = 0; j < YDIM; j++)
        {

            for (int i = 0; i < XDIM; i++)
            {
                int index = (k * YDIM + j) * XDIM + i;
                x[index] = (double)i;
                y[index] = (double)j;
                z[index] = (double)k;
                rho[index] = 1;
                ux[index] = U0 * sin(i * (2 * pi / XDIM)) * (cos(3 * j * (2 * pi / YDIM)) * cos(k * (2 * pi / ZDIM)) - cos(j * (2 * pi / YDIM)) * cos(3 * k * (2 * pi / ZDIM)));
                uy[index] = U0 * sin(j * (2 * pi / YDIM)) * (cos(3 * k * (2 * pi / ZDIM)) * cos(i * (2 * pi / XDIM)) - cos(k * (2 * pi / ZDIM)) * cos(3 * i * (2 * pi / XDIM)));
                uz[index] = U0 * sin(k * (2 * pi / ZDIM)) * (cos(3 * i * (2 * pi / XDIM)) * cos(j * (2 * pi / YDIM)) - cos(i * (2 * pi / XDIM)) * cos(3 * j * (2 * pi / YDIM)));
            }
        }
    }


    //  Create the ordered Zone
    INTEGER4 ZoneType = 0;
    I = TECZNE142((char*)"Ordered Zone",
        &ZoneType,
        &IMax,
        &JMax,
        &KMax,
        &ICellMax,
        &JCellMax,
        &KCellMax,
        &SolTime,
        &StrandID,
        &unused,
        &IsBlock,
        &NFConns,
        &FNMode,
        &TotalNumFaceNodes,
        &TotalNumBndryFaces,
        &TotalNumBndryConnections,
        NULL, // PassiveVarList
        valueLocations, // ValueLocation
        NULL, // ShareVarFromZone
        &ShrConn);

    /*
     * The nodal index ranges of each partition.
     * Interior boundary index ranges must coincide. 
     */
    INTEGER4 dim[3];
    dim[0] = 2; dim[1] = 2; dim[2] = 1;//please specify the division of processors in X,Y and Z respectively
    
    INTEGER4 nx = 0, ny = 0, nz = 0;

    INTEGER4 local_nx = XDIM / dim[0];// X range for a processor
    INTEGER4 local_ny = YDIM / dim[1];// Y range for a processor
    INTEGER4 local_nz = ZDIM / dim[2];// Z range for a processor

    INTEGER4 xmin = 1, ymin = 1, zmin = 1, xmax = 1, ymax = 1, zmax = 1;
    INTEGER4 coord[3] = { 0,0,0 };
    
    
    INTEGER4 numPartitions = dim[0] * dim[1] * dim[2];
    int numProcesses = dim[0] * dim[1] * dim[2];
    INTEGER4 partitionIndices[4][6];//Instead of 4 type the total number of processors

    //The following loop prints the partitions to be printed according to the coordinates
    for (int t = 0; t < numProcesses; t++) {
        if (nx == dim[0]) {
            nx = 0; ny++;
        }
        if (ny == dim[0]) {
            ny = 0; nz++;
        }
        if (nx==0) {
            xmin = 1+ nx * local_nx;
            xmax = xmin + (local_nx-1 );
        }

        else {
            xmin = nx * local_nx;
            xmax = xmin + (local_nx);
        }
        if (ny==0) {
            ymin = 1 + ny * local_ny;
            ymax = ymin + (local_ny - 1);
        }
        else {
            ymin = ny * local_ny;
            ymax = ymin + (local_ny);
        }
        if (nz==0) {
            zmin = 1 + nz * local_nz;
            zmax = zmin + (local_nz - 1);
        }
        else {
            zmin = nz * local_nz;
            zmax = zmin + (local_nz);
        }
        
        partitionIndices[t][0] = xmin;
        partitionIndices[t][1] = ymin;
        partitionIndices[t][2] = zmin;
        partitionIndices[t][3] = xmax;
        partitionIndices[t][4] = ymax;
        partitionIndices[t][5] = zmax;
        
        nx++;
    }
    

// Output partitions
    INTEGER4 *partitionOwners= NULL;
    partitionOwners= (INTEGER4*)malloc(numPartitions *sizeof(INTEGER4));
    for (INTEGER4 ptn = 0; ptn < numPartitions; ++ptn)
        partitionOwners[ptn]=(ptn % commSize);
    TECZNEMAP142(&numPartitions, &partitionOwners[0]);

    for (INTEGER4 partition = 1; partition <= numPartitions; ++partition)
    {

        if (partitionOwners[partition - 1] == commRank || commRank == mainRank)
        {

            INTEGER4 partitionIMin = partitionIndices[partition - 1][0];
            INTEGER4 partitionJMin = partitionIndices[partition - 1][1];
            INTEGER4 partitionKMin = partitionIndices[partition - 1][2];
            INTEGER4 partitionIMax = partitionIndices[partition - 1][3];
            INTEGER4 partitionJMax = partitionIndices[partition - 1][4];
            INTEGER4 partitionKMax = partitionIndices[partition - 1][5];
            I = TECIJKPTN142(&partition, &partitionIMin, &partitionJMin, &partitionKMin, &partitionIMax, &partitionJMax, &partitionKMax);
            I = outputVarData(x, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(y, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(z, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(rho, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(ux, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(uy, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
            I = outputVarData(uz, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
        }

    }



    I = TECEND142();

    MPI_Finalize();


    delete[] x;
    delete[] y;
    delete[] z;
    delete[] rho;
    delete[] ux;
    delete[] uy;
    delete[] uz;

    return 0;
}

INTEGER4 outputVarData(
    double* var,
    INTEGER4 iDim, INTEGER4 jDim, INTEGER4 kDim,
    INTEGER4 iMin, INTEGER4 jMin, INTEGER4 kMin,
    INTEGER4 iMax, INTEGER4 jMax, INTEGER4 kMax)
{
    INTEGER4 count = iMax - iMin + 1;
    INTEGER4 isDouble = 1;
    INTEGER4 result = 0;
    for (INTEGER4 k = kMin; result == 0 && k <= kMax; ++k)
    {
        for (INTEGER4 j = jMin; result == 0 && j <= jMax; ++j)
        {
            INTEGER4 index = ((k - 1) * jDim + j - 1) * iDim + iMin - 1;
            result = TECDAT142(&count, &var[index], &isDouble);
        }
    }
    return result;
}
