// This example creates an IJK-ordered zone in 4 partitions.

#if defined TECIOMPI
#include "mpi.h"
#endif

using namespace std;
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "TECIO.h"
#define pi acos(-1)

#ifndef NULL
#define NULL 0
#endif

#define XDIM 100
#define YDIM 100
#define ZDIM 100

//#define IS_DOUBLE



INTEGER4 outputVarData(
    double* var,
    INTEGER4 iDim, INTEGER4 jDim, INTEGER4 kDim,
    INTEGER4 iMin, INTEGER4 jMin, INTEGER4 kMin,
    INTEGER4 iMax, INTEGER4 jMax, INTEGER4 kMax);

int main(int argc, char** argv)
{
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 1;
    INTEGER4 FileType   = 0;
    INTEGER4 FileFormat = 1; // SZPLT; .PLT not supported for partitioned zones
    INTEGER4 I          = 0; // Used to track return codes
    #if defined TECIOMPI
        MPI_Init(&argc, &argv);
        MPI_Comm mpiComm = MPI_COMM_WORLD;
        int commSize;
        MPI_Comm_size(mpiComm, &commSize);
        int commRank;
        MPI_Comm_rank(mpiComm, &commRank);
        #if defined _DEBUG
            if (commRank == 0)
            {
                cout << "Press return to continue" << endl;
                cin.get();
            }
        #endif
    #endif

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"IJK Ordered Zone",
                  (char*)"X Y Z RHO UX UY UZ",
                  (char*)"x.szplt",
                  (char*)".",
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    #if defined TECIOMPI
        INTEGER4 mainRank = 0;
        I = TECMPIINIT142(&mpiComm, &mainRank);
    #endif

    double* x = new double[XDIM * YDIM * ZDIM];
    double* y = new double[XDIM * YDIM * ZDIM];
    double* z = new double[XDIM * YDIM * ZDIM];
    double* rho = new double[XDIM * YDIM * ZDIM];
    double* ux = new double[XDIM * YDIM * ZDIM];
    double* uy = new double[XDIM * YDIM * ZDIM];
    double* uz = new double[XDIM * YDIM * ZDIM];

    INTEGER4 IMax                     = XDIM;
    INTEGER4 JMax                     = YDIM;
    INTEGER4 KMax                     = ZDIM;
    INTEGER4 ICellMax                 = 0;
    INTEGER4 JCellMax                 = 0;
    INTEGER4 KCellMax                 = 0;
    INTEGER4 DIsDouble                = 1;
    double   SolTime                  = 360.0;
    INTEGER4 StrandID                 = 0;      /* StaticZone */
    INTEGER4 unused                   = 0;      // ParentZone is no longer used
    INTEGER4 IsBlock                  = 1;      /* Block */
    INTEGER4 NFConns                  = 0;
    INTEGER4 FNMode                   = 0;
    INTEGER4 TotalNumFaceNodes        = 1;
    INTEGER4 TotalNumBndryFaces       = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 valueLocations[]         = {1, 1, 1, 1, 1, 1, 1};
    INTEGER4 ShrConn                  = 0;
    /*
    string u;
    double c;
    ifstream file;
    file.open("C:\\Users\\Trishant Pamnani\\Downloads\\Flow_field_20.dat");
    if (!file.is_open())
    {
        cout << "error";
        return 1;
    }
    for (int l = 0; l < 13; l++) {
        file >> u;
    }*/
    double U0 = 0.05;
    /*ux[i][j][k] = U0*sin(i*(2*pi/nx))*(cos(3*j*(2*pi/ny))*cos(k*(2*pi/nz)) - cos(j*(2*pi/ny))*cos(3*k*(2*pi/nz)));

                  uy[i][j][k] = U0*sin(j*(2*pi/ny))*(cos(3*k*(2*pi/nz))*cos(i*(2*pi/nx)) - cos(k*(2*pi/nz))*cos(3*i*(2*pi/nx)));

                  uz[i][j][k] = U0*sin(k*(2*pi/nz))*(cos(3*i*(2*pi/nx))*cos(j*(2*pi/ny)) - cos(i*(2*pi/nx))*cos(3*j*(2*pi/ny)));
              */
    for (int k = 0; k < ZDIM; k++)
    {

        for (int j = 0; j < YDIM; j++)
        {

            for (int i = 0; i < XDIM; i++)
            {   
                int index = (k * YDIM + j) * XDIM + i;
                x[index] = (double)i ;
                y[index] = (double)j ;
                z[index] = (double)k ;
                //file >> c;
                //file >> c;
                //file >> c;

                //file >> c;
                rho[index] = 1;
                //file >> c;
                ux[index] = U0 * sin(i * (2 *  pi / XDIM)) * (cos(3 * j * (2 * pi / YDIM)) * cos(k * (2 * pi / ZDIM)) - cos(j * (2 * pi / YDIM)) * cos(3 * k * (2 * pi / ZDIM)));
                //file >> c;
                uy[index] = U0 * sin(j * (2 * pi / YDIM)) * (cos(3 * k * (2 * pi / ZDIM)) * cos(i * (2 * pi / XDIM)) - cos(k * (2 * pi / ZDIM)) * cos(3 * i * (2 * pi / XDIM)));
                //file >> c;
                uz[index] = U0 * sin(k * (2 * pi / ZDIM)) * (cos(3 * i * (2 * pi / XDIM)) * cos(j * (2 * pi / YDIM)) - cos(i * (2 * pi / XDIM)) * cos(3 * j * (2 * pi / YDIM)));
            }
        }
    }
/*
    // Ordered Zone Parameters
    for(int i = 0; i < XDIM; ++i)
        for(int j = 0; j < YDIM; ++j)
            for(int k = 0; k < ZDIM; ++k)
            {
                int index = (k * YDIM + j) * XDIM + i;
                x[index] = (float)i;
                y[index] = (float)j;
                z[index] = (float)k;
                if (i < XDIM - 1 && j < YDIM - 1 && k < ZDIM - 1)
                {
                    int cindex = (k * (YDIM - 1) + j) * (XDIM - 1) + i;
                    p[cindex] = (x[index] + 0.5f) * (y[index] + 0.5f) * (z[index] + 0.5f);
                }
            }
*/
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
     * Interior boundary index ranges must coincide. For example, if one partition's I index
     * range goes from 1 to 3, its neighboring partition must start with an I index of 3.
     */
    INTEGER4 partitionIndices[4][6] = {
        { 1, 1, 1, XDIM / 3 + 1, YDIM, ZDIM },
        { XDIM / 3 + 1, 1, 1, XDIM, YDIM / 3 + 1, ZDIM },
        { XDIM / 3 + 1, YDIM / 3 + 1, 1, XDIM, YDIM, ZDIM / 2 + 1 },
        { XDIM / 3 + 1, YDIM / 3 + 1, ZDIM / 2 + 1, XDIM, YDIM, ZDIM }
    };

    // Output partitions
    #if defined TECIOMPI
        INTEGER4 numPartitions = 4;
        std::vector<INTEGER4> partitionOwners;
        for (INTEGER4 ptn = 0; ptn < numPartitions; ++ptn)
            partitionOwners.push_back(ptn % commSize);
        TECZNEMAP142(&numPartitions, &partitionOwners[0]);
    #endif

    for (INTEGER4 partition = 1; partition <= 4; ++partition)
    {
        #if defined TECIOMPI
            if (commRank == mainRank || partitionOwners[partition - 1] == commRank)
            {
        #endif
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

        #if defined TECIOMPI
            }
        #endif
    }

    // The I=IMax boundary is output as three unpartitioned surface zones.
    // For MPI, each of these is output by the MPI rank that also output
    // the partition that it boundaries (partitions 2 through 4 touch the
    // I=IMax boundary).
    /*for (INTEGER4 partition = 1; partition <= 4; ++partition)
    {
        #if defined TECIOMPI
            // Only the main rank and the rank that owns the zone need call TECZNE.
            // The others may, but it becomes a no-op for them.
            if (commRank == mainRank || partitionOwners[partition - 1] == commRank)
            {
        #endif
        INTEGER4 surfaceZoneJMin = partitionIndices[partition - 1][1];
        INTEGER4 surfaceZoneKMin = partitionIndices[partition - 1][2];
        INTEGER4 surfaceZoneJMax = partitionIndices[partition - 1][4];
        INTEGER4 surfaceZoneKMax = partitionIndices[partition - 1][5];
        IMax = 1;
        JMax = surfaceZoneJMax - surfaceZoneJMin + 1;
        KMax = surfaceZoneKMax - surfaceZoneKMin + 1;
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
        #if defined TECIOMPI
            INTEGER4 one = 1;
            I = TECZNEMAP142(&one, &partitionOwners[partition - 1]);
        #endif
        I = outputVarData(x, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        I = outputVarData(y, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        I = outputVarData(z, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        // Output the cell-centered values for the iMax - .5 plane of cells
        I = outputVarData(p, XDIM - 1, YDIM - 1, ZDIM - 1,
            XDIM - 1, surfaceZoneJMin, surfaceZoneKMin, XDIM - 1, surfaceZoneJMax - 1, surfaceZoneKMax - 1);
        #if defined TECIOMPI
            }
        #endif
    }*/

    I = TECEND142();
    #if defined TECIOMPI
        MPI_Finalize();
    #endif

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
            INTEGER4 index = ((k-1) * jDim + j -1) * iDim + iMin-1;
            result = TECDAT142(&count, &var[index], &isDouble);
        }
    }
    return result;
}

