#include "epot_bicgstabsolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "geomplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include <cmath>
#include <string>
#include <iostream>
#include <readascii.hpp>

double U = 25.0e3; //не больше 30 кВ, а лучше 25кВ
double ground = -40.0e3;
double  const cond_length = 0.15;
double const collector_length = 0.1;

bool cyl_1(double x,double y,double z){
        return (x*x + y*y >= (0.04)*(0.04) && z <= 0.8);
    }
bool cond_down(double x,double y,double z){
        return(z>=0.82 && z<= (0.82+ cond_length ) && y>=-0.07 && y<=0.07 && x>=-0.070 && x <= -0.065); 
    }
bool cond_up(double x,double y,double z){
        return(z>=0.82 && z<=(0.82 + cond_length ) && y>=-0.07 && y<=0.07 && x >= 0.065 && x<=0.070); 
    }
bool cyl_2(double x,double y,double z){
        return (x*x + y*y >= (0.065)*(0.065) && z >= (0.84+ cond_length) && z<=(0.84+ cond_length + collector_length)); 
    }
bool cyl_3(double x,double y,double z){
        return (x*x + y*y >= (0.04)*(0.04) && z >= (0.84 + cond_length + collector_length) && z<=1.2); 
    }



void simu(void){

    double const x = 140;
    double const y = 140;
    double const z = 450;
    double const mesh_step = 0.002;
    Geometry geom( MODE_3D, Int3D(ceil(x/(mesh_step*1000)),ceil(y/(mesh_step*1000)),ceil(z/(mesh_step*1000))), Vec3D(-0.07,-0.07,0.7), mesh_step );

    
    Solid *Cyl_1 = new FuncSolid( cyl_1 );
    geom.set_solid( 7, Cyl_1 );
    Solid *Cond_down = new FuncSolid( cond_down );
    geom.set_solid( 8, Cond_down );
    Solid *Cond_up = new FuncSolid( cond_up);
    geom.set_solid( 9, Cond_up );
    Solid *Cyl_2 = new FuncSolid( cyl_2 );
    geom.set_solid( 10, Cyl_2 );
    Solid *Cyl_3 = new FuncSolid( cyl_3 );
    geom.set_solid( 11, Cyl_3 );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  ground) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, ground) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, ground-U) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, ground) );
    geom.set_boundary( 11, Bound(BOUND_DIRICHLET, ground) );

    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                    FIELD_EXTRAPOLATE,FIELD_EXTRAPOLATE,
                                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );
    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    ReadAscii din( "particles", 10 );
    GeomPlotter geomplotter( geom );
    std::string buffer = "plot";

    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics;
    
  

    for( int i = 0; i < 5; i++ ) {
        solver.solve( epot, scharge );
        efield.recalculate();
        pdb.clear();
        buffer = "plot" ;
        for( size_t j = 0; j < din.rows(); j+=1) {
        double I  = din[0][j];
        double m  = din[1][j];
        double q  = din[2][j];
        double t = din[3][j];
        double x = din[4][j];
        double vx = din[5][j];
        double y = din[6][j];
        double vy = din[7][j];
        double z = din[8][j];
        double vz = din[9][j];

        pdb.add_particle( I, q/CHARGE_E, m/MASS_U, ParticleP3D(t,x,vx,y,vy,z,vz) );
    }
        buffer = buffer + to_string(i) + ".png";
   pdb.iterate_trajectories( scharge, efield, bfield );
  
    geomplotter.set_size( 750, 750 );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_view(VIEW_ZX);
    geomplotter.plot_png( buffer );
    }
    diagnostics.push_back(DIAG_CURR);    
    pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diagnostics );
    pdb.trajectories_at_plane( tdata, AXIS_Z, 0.82, diagnostics );

}

int main( int argc, char **argv )
{
    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
        simu();
    } catch ( Error e ) {
        e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
