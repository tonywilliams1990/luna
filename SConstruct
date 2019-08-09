import os.path
import os
import glob


# compiler
c_comp = 'g++'
#c_comp = 'mpic++'

# library names
Luna_lib = 'Luna'

# default output format for messaging
red = "\033[1;31m"
yellow = "\033[1;33m"
green = "\033[1;32m"
blue = "\033[1;34m"
off = "\033[0m"

def message( col, text ):
    print col + " * " + text + off

# Output message to user
message( blue, " -----  Building ----- ")

# set the build dir
topdir     = os.getcwd()
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
libs_str   = Luna_lib + ' '
preproc    = ''
#opts = ' -O2 -std=c++17 -Wall -Wextra '
#opts = ' -O2 -std=c++17 -fopenmp '
#opts = ' -Ofast -std=c++17 -fopenmp '
opts = ' -O3 -std=c++17 -funroll-loops -fopenmp '
#opts = ' -O3 -std=c++11'# -fopenmp'
link_flags = ' -fopenmp '

# Set the rpath for the linker to find petsc/slepc-3
rpath = []

# source is all cpp files, so let's glob them
src = glob.glob('src/*.cpp')

# Split the strings into list elements for lib/inc directoris
incdir = incdir_str.split()
libdir = libdir_str.split()

# Initialise the environment
env = Environment( CXX = c_comp, CPPPATH = incdir, CCFLAGS = opts + preproc, LINKFLAGS = link_flags, LIBPATH = libdir)

# Now check the environment is complete
conf = Configure( env )
env = conf.Finish()

libs = libs_str.split()

# Build the libraries in ./lib
env.StaticLibrary('lib/' + Luna_lib, [src] )

SConscript('Examples/SConscript', exports='env opts preproc link_flags rpath incdir libdir topdir libs' )
