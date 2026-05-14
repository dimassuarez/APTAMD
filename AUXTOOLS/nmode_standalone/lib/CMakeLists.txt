#Amber utility library
#In the Makefile build system, parts of this library are built for various modules as needed, but in the CMake system we build all of them now
#and link in the single library to the modules
project(amber_common)

set(LIB_FORTRAN_SOURCES sys.F veclib.F90 random.F90 amopen.F mexit.F
 nxtsec.F rgroup.F rfree.F matinv.F constants.F90 AmberNetcdf.F90 ligenb.F util.F)
 
set(LIB_C_SOURCES wallclock.c binpos.c)

#set optimization flags
apply_optimization_declarations(${LIB_FORTRAN_SOURCES})
set_property(SOURCE wallclock.c PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS})
set_property(SOURCE binpos.c PROPERTY COMPILE_OPTIONS ${NO_OPT_CFLAGS})

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

add_library(amber_common STATIC ${LIB_FORTRAN_SOURCES} ${LIB_C_SOURCES})
target_link_libraries(amber_common netcdff netcdf) 
make_pic_if_needed(amber_common)
config_module_dirs(amber_common ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})

install_libraries(amber_common)


