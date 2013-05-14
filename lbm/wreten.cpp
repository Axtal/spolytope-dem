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

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/dem/domain.h>

struct UserData
{
    Array<Cell *>  xmin0;
    Array<Cell *>  xmax0;
    Array<Cell *>  xmin1;
    Array<Cell *>  xmax1;
    Array<Cell *>  ymin0;
    Array<Cell *>  ymax0;
    Array<Cell *>  ymin1;
    Array<Cell *>  ymax1;
    Array<Cell *>  zmin0;
    Array<Cell *>  zmax0;
    Array<Cell *>  zmin1;
    Array<Cell *>  zmax1;
    size_t         block;       ///< Dimensions of the checkerboard block
    double          Head;       ///< Current hydraulic head
    double          Orig;       ///< Original hydraulic head
    double            Tf;
    double            Kn;
    double           ome;
    double         dtOut;
    double          time;
    double           rho;
    Vec3_t            Dp;       ///< Vector with the Pc gradient
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.Time>dat.time)
    {
        dat.time += dat.dtOut;
    }
    double a   = M_PI/dat.ome;
    double rho = dat.Head*((1.0/a)*(dat.time-a*(floor(dat.time/a)+0.5))*pow(-1.0,floor(dat.time/a))+0.5)+dat.Orig;
    double rho0min;
    double rho1min;
    double rho0max;
    double rho1max;
    
    if (dat.block<2)
    { 
        if (fabs(dat.Dp(0))>1.0e-12)
        {
            rho0min = 0.999*((rho - dat.rho)*dat.Dp(0) + dat.rho);
            rho1min = 0.001*((rho - dat.rho)*dat.Dp(0) + dat.rho);
            rho0max = 0.001*((dat.rho - rho)*dat.Dp(0) + dat.rho);
            rho1max = 0.999*((dat.rho - rho)*dat.Dp(0) + dat.rho);
            for (size_t i=0;i<dat.xmin0.Size();i++)
            {
                Cell * c = dat.xmin0[i];
                c->RhoBC = rho0min;
                c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
                c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.xmin1[i];
                c->RhoBC = rho1min;
                c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
                c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.xmax0[i];
                c->RhoBC = rho0max;
                c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
                c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
                c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);
                
                c = dat.xmax1[i];
                c->RhoBC = rho1max;
                c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
                c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
                c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);
            }
        }
        if (fabs(dat.Dp(1))>1.0e-12)
        {
            rho0min = 0.999*((rho - dat.rho)*dat.Dp(1) + dat.rho);
            rho1min = 0.001*((rho - dat.rho)*dat.Dp(1) + dat.rho);
            rho0max = 0.001*((dat.rho - rho)*dat.Dp(1) + dat.rho);
            rho1max = 0.999*((dat.rho - rho)*dat.Dp(1) + dat.rho);
            for (size_t i=0;i<dat.ymin0.Size();i++)
            {
                Cell * c = dat.ymin0[i];
                c->RhoBC = rho0min;
                c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
                c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.ymin1[i];
                c->RhoBC = rho1min;
                c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
                c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.ymax0[i];
                c->RhoBC = rho0max;
                c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
                c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.ymax1[i];
                c->RhoBC = rho1max;
                c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
                c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);
            }
        }
        if (fabs(dat.Dp(2))>1.0e-12)
        {
            rho0min = 0.999*((rho - dat.rho)*dat.Dp(2) + dat.rho);
            rho1min = 0.001*((rho - dat.rho)*dat.Dp(2) + dat.rho);
            rho0max = 0.001*((dat.rho - rho)*dat.Dp(2) + dat.rho);
            rho1max = 0.999*((dat.rho - rho)*dat.Dp(2) + dat.rho);
            for (size_t i=0;i<dat.zmin0.Size();i++)
            {
                Cell * c = dat.zmin0[i];
                c->RhoBC = rho0min;
                c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
                c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.zmin1[i];
                c->RhoBC = rho1min;
                c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
                c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.zmax0[i];
                c->RhoBC = rho0max;
                c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
                c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);

                c = dat.zmax1[i];
                c->RhoBC = rho1max;
                c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
                c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
                c->Rho = c->VelDen(c->Vel);
            }
        }
    }
    else
    {
        for (size_t i=0;i<dat.xmin0.Size();i++)
        {
            Cell * c = dat.xmin0[i];
            size_t chkboard = ((c->Index(1)/dat.block)%2 + (c->Index(2)/dat.block)%2)%2;
            if (chkboard==0)
            {
                rho0min = 0.999*rho;
                rho1min = 0.001*rho;
                rho0max = 0.001*(2.0*dat.rho - rho);
                rho1max = 0.999*(2.0*dat.rho - rho);
            }
            else if (chkboard==1)
            {
                rho0max = 0.999*rho;
                rho1max = 0.001*rho;
                rho0min = 0.001*(2.0*dat.rho - rho);
                rho1min = 0.999*(2.0*dat.rho - rho);
            }

            c->RhoBC = rho0min;
            //c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            //c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.xmin1[i];
            c->RhoBC = rho1min;
            //c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            //c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.xmax0[i];
            c->RhoBC = rho0max;
            //c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            //c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            //c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);
            
            c = dat.xmax1[i];
            c->RhoBC = rho1max;
            //c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            //c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            //c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);
        }

        for (size_t i=0;i<dat.ymin0.Size();i++)
        {
            Cell * c = dat.ymin0[i];
            size_t chkboard = ((c->Index(0)/dat.block)%2 + (c->Index(2)/dat.block)%2)%2;
            if (chkboard==0)
            {
                rho0min = 0.999*rho;
                rho1min = 0.001*rho;
                rho0max = 0.001*(2.0*dat.rho - rho);
                rho1max = 0.999*(2.0*dat.rho - rho);
            }
            else if (chkboard==1)
            {
                rho0max = 0.999*rho;
                rho1max = 0.001*rho;
                rho0min = 0.001*(2.0*dat.rho - rho);
                rho1min = 0.999*(2.0*dat.rho - rho);
            }

            c->RhoBC = rho0min;
            //c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
            //c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.ymin1[i];
            c->RhoBC = rho1min;
            //c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
            //c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.ymax0[i];
            c->RhoBC = rho0max;
            //c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
            //c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.ymax1[i];
            c->RhoBC = rho1max;
            //c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
            //c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);
        }
    
        for (size_t i=0;i<dat.zmin0.Size();i++)
        {
            Cell * c = dat.zmin0[i];
            size_t chkboard = ((c->Index(0)/dat.block)%2 + (c->Index(1)/dat.block)%2)%2;
            if (chkboard==0)
            {
                rho0min = 0.999*rho;
                rho1min = 0.001*rho;
                rho0max = 0.001*(2.0*dat.rho - rho);
                rho1max = 0.999*(2.0*dat.rho - rho);
            }
            else if (chkboard==1)
            {
                rho0max = 0.999*rho;
                rho1max = 0.001*rho;
                rho0min = 0.001*(2.0*dat.rho - rho);
                rho1min = 0.999*(2.0*dat.rho - rho);
            }

            c->RhoBC = rho0min;
            //c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
            //c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.zmin1[i];
            c->RhoBC = rho1min;
            //c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
            //c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.zmax0[i];
            c->RhoBC = rho0max;
            //c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
            //c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);

            c = dat.zmax1[i];
            c->RhoBC = rho1max;
            //c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
            //c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
            //c->Rho = c->VelDen(c->Vel);
            c->Initialize(c->RhoBC,OrthoSys::O);
        }
    }
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double water = 0.0;
    double oil   = 0.0;
    double Sr    = 0.0;
    size_t nw    = 0;
    size_t no    = 0;
    for (size_t i=0;i<dom.Lat[1].Cells.Size();i++)
    {
        double wr = dom.Lat[1].Cells[i]->Rho;
        double ar = dom.Lat[0].Cells[i]->Rho;
        if (dom.Lat[1].Cells[i]->IsSolid) continue;
        if (wr>0.5*dat.rho) 
        {
            Sr+=1.0;
            water+=(wr + ar + dom.Gmix*wr*ar)/3.0;
            nw++;
        }
        if (ar>0.5*dat.rho)
        {
            oil  +=(wr + ar + dom.Gmix*wr*ar)/3.0;
            no++;
        }
    }
    Sr/=(dom.Lat[0].Cells.Size()*(1-dom.Lat[0].SolidFraction()));
    if (nw>0) water/=nw;
    if (no>0) oil  /=no;
    double rhow = 0.0;
    double rhoo = 0.0;
    size_t nfb  = 0;
    size_t nfo  = 0;
    for (size_t i=0;i<dom.Lat[0].Ndim(1);i++)
    for (size_t j=0;j<dom.Lat[0].Ndim(2);j++)
    {
        Cell * c = dom.Lat[0].GetCell(iVec3_t(1,i,j));
        if (c->IsSolid) continue;
        rhow += c->Rho;        
        nfb++;
        c = dom.Lat[1].GetCell(iVec3_t(dom.Lat[1].Ndim(0)-2,i,j));
        rhoo += c->Rho;        
        nfo++;
    }
    rhow/=nfb;
    rhoo/=nfo;
    double Pc; 
    double rho;
    
    double a   = M_PI/dat.ome;
    rho        = dat.Head*((1.0/a)*(dat.time-a*(floor(dat.time/a)+0.5))*pow(-1.0,floor(dat.time/a))+0.5)+dat.Orig;
    Pc         = (2.0*(rho - dat.rho) + dom.Gmix*(rho*rho*0.999*0.001 - (2.0*dat.rho - rho)*(2.0*dat.rho - rho)*0.999*0.001))/3.0;
    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << rhoo << Util::_8s << rhow << Util::_8s << water << Util::_8s << oil << Util::_8s << Pc << Util::_8s << Sr << std::endl;
}

int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc==3) Nproc = atoi(argv[2]);

    String fileDEM;
    String fileLBM;
    bool   Render   = true;
    size_t N        = 200;
    double Gs0      = -0.53;
    double Gs1      = -0.53;
    double Gmix     = 2.0;
    double nu       = 0.05;
    double dt       = 1.0;
    double Tf       = 10000.0;
    double dtOut    = 50.0;
    double HeadStep = 1000.0;
    double rho      = 200.0;
    double ome      = 2.0;
    double Head     = 500.0;
    double Orig     = 54.0;
    size_t oct      = 1;
    double DPx      = 1.0;
    double DPy      = 1.0;
    double DPz      = 1.0;
    int    outlimit = 1;
    size_t buffer   = 1;
    {
        infile >> fileDEM;   infile.ignore(200,'\n');
        infile >> fileLBM;   infile.ignore(200,'\n');
        infile >> Render;    infile.ignore(200,'\n');
        infile >> N;         infile.ignore(200,'\n');
        infile >> Gs0;       infile.ignore(200,'\n');
        infile >> Gs1;       infile.ignore(200,'\n');
        infile >> Gmix;      infile.ignore(200,'\n');
        infile >> nu;        infile.ignore(200,'\n');
        infile >> dt;        infile.ignore(200,'\n');
        infile >> Tf;        infile.ignore(200,'\n');
        infile >> dtOut;     infile.ignore(200,'\n');
        infile >> HeadStep;  infile.ignore(200,'\n');
        infile >> rho;       infile.ignore(200,'\n'); 
        infile >> ome;       infile.ignore(200,'\n');
        infile >> Head;      infile.ignore(200,'\n');
        infile >> Orig;      infile.ignore(200,'\n');
        infile >> oct;       infile.ignore(200,'\n'); 
        infile >> DPx;       infile.ignore(200,'\n');
        infile >> DPy;       infile.ignore(200,'\n');
        infile >> DPz;       infile.ignore(200,'\n');
        infile >> outlimit;  infile.ignore(200,'\n'); 
        infile >> buffer;    infile.ignore(200,'\n'); 
    }
    Array<double> nua(2);
    nua[0] = nu;
    nua[1] = nu;


    DEM::Domain DemDom;
    DemDom.Load(fileDEM.CStr());
    Array<int> idx(6);
    idx = -2,-3,-4,-5,-6,-7;
    DemDom.DelParticles(idx);
    Vec3_t Xmin,Xmax;
    DemDom.BoundingBox(Xmin,Xmax);
    //int    bound = -2;
    int    bound = outlimit;
    double dx = (Xmax(0)-Xmin(0))/(N-2*bound);
    size_t Ny = (Xmax(1)-Xmin(1))/dx + 2*bound;
    size_t Nz = (Xmax(2)-Xmin(2))/dx + 2*bound;
    DemDom.Center(0.5*(Xmax-Xmin)+Vec3_t(bound*dx,bound*dx,bound*dx));
    LBM::Domain Dom(D3Q15, nua, iVec3_t(N,Ny,Nz), 1.0, 1.0);
    Dom.PrtVec = false;
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<Ny;j++)
        {
            for (int k=0;k<Nz;k++)
            {
                Vec3_t pos((i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx);
                for (size_t n=0;n<DemDom.Particles.Size();n++)
                {
                    DEM::Particle *P = DemDom.Particles[n];
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                }
            }
        }
    }

    UserData dat;
    Dom.UserData = &dat;
    if (oct==0)
    {
        double Giso = DPx;
        double Gdev = DPy;
        double th   = DPz;
        DPx    = Giso/sqrt(3.0) + sqrt(2.0/3.0)*Gdev*sin(M_PI*th/180.0-2.0*M_PI/3.0);
        DPy    = Giso/sqrt(3.0) + sqrt(2.0/3.0)*Gdev*sin(M_PI*th/180.0);
        DPz    = Giso/sqrt(3.0) + sqrt(2.0/3.0)*Gdev*sin(M_PI*th/180.0+2.0*M_PI/3.0);
    }

    dat.Tf       = Tf;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.Orig     = Orig;
    dat.dtOut    = HeadStep;
    dat.time     = 0.0;
    dat.rho      = rho;
    dat.Head     = Head;
    dat.Dp       = Vec3_t(DPx,DPy,DPz);
    dat.Dp      /= norm(dat.Dp);
    dat.block    = oct;

    Dom.Lat[0].G = 0.0;
    Dom.Lat[0].Gs= Gs0;
    Dom.Lat[1].G = 0.0;
    Dom.Lat[1].Gs= Gs1;
    Dom.Gmix     = Gmix;

	// set inner obstacles
    
    for (int i=0;i<Ny;i++)
    for (int j=0;j<Nz;j++)
    {
        dat.xmin0.Push(Dom.Lat[0].GetCell(iVec3_t(0  ,i,j)));
        dat.xmax0.Push(Dom.Lat[0].GetCell(iVec3_t(N-1,i,j)));
        dat.xmin1.Push(Dom.Lat[1].GetCell(iVec3_t(0  ,i,j)));
        dat.xmax1.Push(Dom.Lat[1].GetCell(iVec3_t(N-1,i,j)));
    }
    for (int i=0;i<N ;i++)
    for (int j=0;j<Nz;j++)
    {
        dat.ymin0.Push(Dom.Lat[0].GetCell(iVec3_t(i,0   ,j)));
        dat.ymax0.Push(Dom.Lat[0].GetCell(iVec3_t(i,Ny-1,j)));
        dat.ymin1.Push(Dom.Lat[1].GetCell(iVec3_t(i,0   ,j)));
        dat.ymax1.Push(Dom.Lat[1].GetCell(iVec3_t(i,Ny-1,j)));
    }
    for (int i=0;i<N ;i++)
    for (int j=0;j<Ny;j++)
    {
        dat.zmin0.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0   )));
        dat.zmax0.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,Nz-1)));
        dat.zmin1.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,0   )));
        dat.zmax1.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,Nz-1)));
    }
    for (int i=0;i<N ;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,Ny-1,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,Nz-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,Ny-1,Nz-1))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0   ))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,Ny-1,0   ))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,Nz-1))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,Ny-1,Nz-1))->IsSolid = true;
    }
    for (int i=0;i<Ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0  ,i,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,i,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,i,Nz-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,i,Nz-1))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0  ,i,0   ))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(N-1,i,0   ))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0  ,i,Nz-1))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(N-1,i,Nz-1))->IsSolid = true;
    }
    for (int i=0;i<Nz;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0  ,0   ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,0   ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,Ny-1,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,Ny-1,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0  ,0   ,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(N-1,0   ,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0  ,Ny-1,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(N-1,Ny-1,i))->IsSolid = true;
    }

    bound = buffer;
    //Initializing values
    if (Util::FileExists(fileLBM))
    {
        hid_t file_id;
        file_id = H5Fopen(fileLBM.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
        float * Density0   = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Vvec0      = new float[3*Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Density1   = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Vvec1      = new float[3*Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        H5LTread_dataset_float(file_id,"Density_0" ,Density0);
        H5LTread_dataset_float(file_id,"Velocity_0",Vvec0   );
        H5LTread_dataset_float(file_id,"Density_1" ,Density1);
        H5LTread_dataset_float(file_id,"Velocity_1",Vvec1   );
        for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
        {
            Cell * c = Dom.Lat[0].Cells[i];
            Vec3_t V;
            V(0) =  Vvec0[3*i  ];
            V(1) =  Vvec0[3*i+1];
            V(2) =  Vvec0[3*i+2];
            c->Initialize(Density0[i],V);
            c = Dom.Lat[1].Cells[i];
            V(0) =  Vvec1[3*i  ];
            V(1) =  Vvec1[3*i+1];
            V(2) =  Vvec1[3*i+2];
            c->Initialize(Density1[i],V);
        }
        delete [] Density0;
        delete [] Vvec0;
        delete [] Density1;
        delete [] Vvec1;
    }
    else
    {
        for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
        {
            Cell * c0 = Dom.Lat[0].Cells[i];
            Cell * c1 = Dom.Lat[1].Cells[i];
            c1->Initialize(0.999*rho, OrthoSys::O);
            c0->Initialize(0.001*rho, OrthoSys::O);
            if (oct<2)
            {
                if ((dat.Dp(0)>1.0e-12)&&(c0->Index(0)<N/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                if ((dat.Dp(1)>1.0e-12)&&(c0->Index(1)<Ny/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                if ((dat.Dp(2)>1.0e-12)&&(c0->Index(2)<Nz/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
            }
            else
            {
                size_t chkboard = ((c0->Index(1)/dat.block)%2 + (c0->Index(2)/dat.block)%2)%2;
                if ((chkboard==0)&&(c0->Index(0)<N/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                else if ((chkboard==1)&&(c0->Index(0)>=(bound-1)*N/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                chkboard = ((c0->Index(0)/dat.block)%2 + (c0->Index(2)/dat.block)%2)%2;
                if ((chkboard==0)&&(c0->Index(1)<Ny/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                else if ((chkboard==1)&&(c0->Index(1)>=(bound-1)*Ny/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                chkboard = ((c0->Index(0)/dat.block)%2 + (c0->Index(1)/dat.block)%2)%2;
                if ((chkboard==0)&&(c0->Index(2)<Nz/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
                else if ((chkboard==1)&&(c0->Index(2)>=(bound-1)*Nz/bound))
                {
                    c0->Initialize(0.999*rho, OrthoSys::O);
                    c1->Initialize(0.001*rho, OrthoSys::O);
                }
            }

            //if      (c0->Index(0)<=abs(bound))
            //{
                //c0->Initialize(0.999*rho, OrthoSys::O);
                //c1->Initialize(0.001*rho, OrthoSys::O);
            //}
            //else if (c0->Index(1)<=abs(bound))
            //{
                //c0->Initialize(0.999*rho, OrthoSys::O);
                //c1->Initialize(0.001*rho, OrthoSys::O);
            //}
            //else if (c0->Index(2)<=abs(bound))
            //{
                //c0->Initialize(0.999*rho, OrthoSys::O);
                //c1->Initialize(0.001*rho, OrthoSys::O);
            //}
            //else
            //{
                //c1->Initialize(0.999*rho, OrthoSys::O);
                //c0->Initialize(0.001*rho, OrthoSys::O);
            //}
        }
    }

    //Solving
    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(),std::ios::out);
    dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "PDen" << Util::_8s << "Rhow" << Util::_8s << "Rhoo" << Util::_8s << "Water" << Util::_8s << "Oil" << Util::_8s << "Pc" << Util::_8s << "Sr" << std::endl;
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH


