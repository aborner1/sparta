###################################
# Trajectory inputs
###################################
variable            V equal 6000.0
variable 	    rho equal 1.13e-6
variable            temp equal 195.5
variable            surftemp equal 1000.0

###################################
# Gas parameters (C02)
###################################
variable            mue equal 1.380E-5
variable            mass equal 73.1E-27
variable            gamma equal 1.28
variable            nden equal ${rho}/${mass}

###################################
# Constants
###################################
variable            boltz equal 1.380658E-23
variable            pi equal 3.14159

###################################
# Simulation initialization standards
###################################
variable            ppc equal 1000
variable            cpmfp equal 1

###################################
# Parameter calculations
###################################
variable            cbar equal sqrt(8./${pi}*${boltz}*${temp}/${mass})
variable            mfp equal 2*${mue}/(${nden}*${mass}*${cbar})
variable            uspeed equal sqrt(${boltz}*${temp}/${mass})

variable            xmin equal 0
variable            xmax equal 3.2
variable            ymin equal 0
variable            ymax equal 6

variable            xncells equal floor((${xmax}-${xmin})/${mfp}*${cpmfp})
variable            yncells equal floor((${ymax}-${ymin})/${mfp}*${cpmfp})

variable            Fnum equal  ${nden}*0.5*(${xmax}-${xmin})*(${ymax}+${ymin})/(${ppc}*0.001)/${xncells}/${yncells}

variable 	    partRatio equal ${nden}/${Fnum}

variable            tstep equal (-${xmin}+${xmax})/${V}/${xncells}/20

###################################
# Print variable values to log file
###################################
print               " Velocity  = ${V}"
print               " Density  = ${rho}"
print               " Number density  = ${nden}"
print               " Temp  = ${temp}"
print               " cbar  = ${cbar}"
print               " mean free path  = ${mfp}"
print               " cells per free stream mean free path = ${cpmfp}"
print               " sound speed  = ${uspeed}"
print               " x-min = ${xmin}"
print               " x-max = ${xmax}"
print               " y-min = 0.0"
print               " y-max = ${ymax}"
print               " x-cells = ${xncells}"
print               " y-cells = ${yncells}"
print               " Timestep = ${tstep}"
print 		    " Fnum = ${Fnum}"
print 		    " Nrho = ${nden}"
print		    " Nrho/Fnum = ${partRatio}"

###################################
# Simulation parameters
###################################
seed	    	    847384
dimension   	    2
global		    nrho ${nden} fnum ${Fnum} gridcut 1e-6 comm/style all
timestep            ${tstep}

###################################
# Grid generation
###################################
boundary	    o ao p
create_box          ${xmin} ${xmax} 0.0 ${ymax} -0.5 0.5
create_grid 	    ${xncells} ${yncells} 1 block * * * 
global              weight cell radius
balance_grid        rcb cell

#####################################
# Gas/Collision Model Specification #
#####################################
species             mars.species N2 CO2 N NO O2 O CO C2 CN C O2+ O+ N+ NO+ CO+ C+ e
mixture             species copy noelectron
mixture             noelectron delete e

mixture		    O2p O2+
mixture		    Op O+
mixture 	    Np N+
mixture 	    NOp NO+
mixture		    COp CO+
mixture 	    Cp C+

mixture             noelectron vstream ${V} 0.0 0.0 temp ${temp}
mixture             noelectron CO2 frac 0.96
mixture             noelectron N2  frac 0.04

#####################################################
# Surface generation and collision specification
#####################################################
read_surf           D10-axi.surf group 1

surf_collide        1 diffuse ${surftemp} 1.0
surf_react          1 prob catalytic.surf
surf_modify         1 collide 1 react 1

fix                 ambi ambipolar e CO+ O+ N+ NO+ O2+ C+

collide             vss species mars.vss relax variable
collide_modify      vremax 20000 yes vibrate smooth rotate smooth

react               tce mars.tce

###################################
# Boundary conditions
###################################
fix                 in emit/face noelectron xlo 

###################################
# Initialize simulation
###################################
create_particles    noelectron n 0
balance_grid        rcb part

###################################
# Unsteady Output
###################################
compute             2 grid all all nrho
compute             4 thermal/grid all all temp
fix                 2 ave/grid all 1 250 250 c_2[*] c_4[*] ave one

compute             1b lambda/grid f_2[1] f_2[2] CO kall

fix                 10 adapt 250 all refine coarsen value c_1b[2] 0.5 2.0 &
                    combine min thresh less more cells 2 2 1
fix                 load balance 250 1.1 rcb time

stats               250
stats_style         step cpu np nattempt ncoll nreact nscoll ngrid maxlevel
stats_modify        flush yes

run                 20000

write_grid          grid.out

collide_modify      vremax 10000 yes
unfix               load
fix                 load balance 10000 1.1 rcb time
unfix               10
unfix               2

###################################
# Steady state Output
###################################
compute             1 grid all COp n nrho tvib
fix                 1 ave/grid all 10 1000 10000 c_1[*] ave running
dump                1 grid all 10000 std_species.* id vol f_1[*]

# Surface pressure and heat flux
compute             5 surf all all press etot
fix                 5 ave/surf all 10 1000 10000 c_5[*] ave running
dump                5 surf all 10000 std_surf.* id f_5[*]

run                 50000
