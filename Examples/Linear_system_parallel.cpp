/// \file  Linear_system.cpp
/// \ingroup Examples
/// Solve some linear systems of equations

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main(int argc, char ** argv)
{
  int mynode, totalnodes;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  if ( mynode == 0 ) // Run this part on a single core
  {

    cout << "----- Linear system parallel -----" << endl;


  }
  // This part is run on every core

  MPI_Barrier( MPI_COMM_WORLD ); // Try to make processes wait but this is not guaranteed

  cout << "Hello world from processor " << mynode << " of " << totalnodes << endl;

  // Delay all processes except 0
  if ( mynode != 0 )
  {
    std::this_thread::sleep_until( std::chrono::system_clock::now()
                                 + std::chrono::seconds(1));
  }

  if ( mynode == 0 )
  {
    cout << "--- FINISHED ---" << endl;
  }
  MPI_Finalize();

}
