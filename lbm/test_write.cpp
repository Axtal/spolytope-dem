/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Torres                                     *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
 
#include <hdf5.h>
#include <hdf5_hl.h>
 
// The number of cells in the X, Y dimensions
#define NX 30
#define NY 30
 
void write_hdf5_data_01()
{
    hid_t     file_id;
    file_id = H5Fcreate("xdmf2d_co.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
 
    // Create the coordinate data.
    double *x = (double *) malloc((NX+1)*(NY+1) * sizeof(double));
    double *y = (double *) malloc((NX+1)*(NY+1) * sizeof(double));
    int ndx = 0;
    for (int j = 0; j < NY+1; j++)
    {
        double yt = double(j) / double(NY);
        double angle = yt * M_PI;
        for (int i = 0; i < NX+1; i++)
        {
            double xt = double(i) / double(NX);
            double R = (1.-xt)*2. + xt*5.;
 
            x[ndx] = R * cos(angle);
            y[ndx] = R * sin(angle);
            ndx++;
        }
    }
 
    hsize_t   dims[1];
    dims[0] = (NY + 1)*(NX + 1);

    H5LTmake_dataset_double(file_id,"X",1,dims,x);
    H5LTmake_dataset_double(file_id,"Y",1,dims,y);
    free(x);
    free(y);
    H5Fclose(file_id);

    file_id = H5Fcreate("xdmf2d_01.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Create the scalar data.
    double *pressure = (double *) malloc(NX*NY * sizeof(double));
 
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int ndx = j * NX + i;
            pressure[ndx] = (double) j;
        }
    }
 
    double *velocityx = (double *) malloc((NX+1)*(NY+1) * sizeof(double));
 
    for (int j = 0; j < NY+1; j++)
    {
        for (int i = 0; i < NX+1; i++)
        {
            int ndx = j * (NX+1) + i;
            velocityx[ndx] = (double) i;
        }
    }
 
    // Write the data file.
    H5LTmake_dataset_double(file_id,"VelocityX",1,dims,velocityx);

    dims[0] = NY*NX;
    H5LTmake_dataset_double(file_id,"Pressure",1,dims,pressure);
 
    free(pressure);
    free(velocityx);
 
}
 
void write_hdf5_data_02()
{
    hid_t     file_id;
    file_id = H5Fcreate("xdmf2d_02.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Create the scalar data.
    double *pressure = (double *) malloc(NX*NY * sizeof(double));
 
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int ndx = j * NX + i;
            pressure[ndx] = 2.0*((double) j);
        }
    }
 
    double *velocityx = (double *) malloc((NX+1)*(NY+1) * sizeof(double));
 
    for (int j = 0; j < NY+1; j++)
    {
        for (int i = 0; i < NX+1; i++)
        {
            int ndx = j * (NX+1) + i;
            velocityx[ndx] = 3.0*((double) i);
        }
    }
    // Write the data file.
    hsize_t   dims[1];
    dims[0] = (NY + 1)*(NX + 1);

    H5LTmake_dataset_double(file_id,"VelocityX",1,dims,velocityx);

    dims[0] = NY*NX;
    H5LTmake_dataset_double(file_id,"Pressure",1,dims,pressure);
 
    free(pressure);
    free(velocityx);
 
    H5Fclose(file_id);
}

void write_xdmf_xml_01()
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen("xdmf2d_01.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d_co.h5:/X\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d_co.h5:/Y\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
    fprintf(xmf, "        xdmf2d_01.h5:/Pressure\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1);
    fprintf(xmf, "        xdmf2d_01.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);

    xmf = fopen("xdmf2d_02.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d_co.h5:/X\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d_co.h5:/Y\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
    fprintf(xmf, "        xdmf2d_02.h5:/Pressure\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1);
    fprintf(xmf, "        xdmf2d_02.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

void write_xdmf_xml_02()
{
    FILE *xmf = 0;
    xmf = fopen("xdmf2d_02.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"%d %d\"/>\n", NY+1, NX+1);
    fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n");
    fprintf(xmf, "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n", (NY+1), (NX+1));
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 2.0 1.0\n", (NY+1), (NX+1));
    fprintf(xmf, "       </DataItem>\n");
    //fprintf(xmf, "     <Geometry GeometryType=\"XY\">\n");
    //fprintf(xmf, "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2 4\"> \n0.0 0.0\n0.0 1.0\n1.0 1.0\n1.0 0.0\n");
    //fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
    fprintf(xmf, "        xdmf2d_02.h5:/Pressure\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1);
    fprintf(xmf, "        xdmf2d_02.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);

}

void write_pvd ()
{
    std::ofstream of("main.pvd", std::ios::out);
    of << "<?xml version=\"1.0\" ?>" << std::endl;
    of << "<XDMFFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    of << "  <Collection>" << std::endl;

    of << "    <DataSet timestep=\"" << 0 << "\" file=\"" << "xdmf2d_01.xmf" << "\" />" << std::endl;
    of << "    <DataSet timestep=\"" << 1 << "\" file=\"" << "xdmf2d_02.xmf" << "\" />" << std::endl;

    of << "  </Collection>" << std::endl;
    of << "</XDMFFile>" << std::endl;
    of.close ();


}

 
int main(int argc, char *argv[])
{
    write_hdf5_data_01();
    write_xdmf_xml_01();
    write_hdf5_data_02();
    write_xdmf_xml_02();
    write_pvd();
 
    return 0;
}
