#include "rabi.h"
#include <fstream>

void progressbar( int percent){
  static int last=-1;
  if(percent==last)
    return;
  last=percent;
  std::string bar;
  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  
  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%  " << std::flush;
}


int main(){
  // Declaration of the objects
  rabi R0,Rd;
  // de-tuning
  double del;

  // delta time for the prints
  double dt=0.01;
  // Final time
  double tf=120;


  // Tuned Rabi system
  R0.init(10,10,0.1);

  // Setting the errors
  R0.Set("rel_error",1e-5);
  R0.Set("abs_error",1e-5);

  cout << "Rabi system with frequency of 10 initialized." << endl;
  cout << "give the value for the detuning: " << endl;
  cin >> del;

  // un-tuned Rabi system
  Rd.init(10,10+del,0.1);
  // Setting the errors
  Rd.Set("rel_error",1e-5);
  Rd.Set("abs_error",1e-5);


  cout << "Computing rabi" << endl;
  ofstream file("rabi.dat");

  // Evolve and save the evolution
  for(double t=0;t<tf;t+=dt){
    progressbar(100*t/tf);
    R0.EvolveSUN(t,t+dt);
    file << t << "\t" << R0.GetExpectationValue(R0.d0,0,0) << "  " 
	 << R0.GetExpectationValue(R0.evol_b0_proj[0],0,0) << "  " 
	 << R0.GetExpectationValue(R0.evol_b0_proj[1],0,0) << endl;
    R0.set_evol();
  }
  file.close();
  file.open("rabi_detuned.dat");
  cout << endl << "Computing detuned rabi" << endl;
  for(double t=0;t<tf;t+=dt){
    progressbar(100*t/tf);
    Rd.EvolveSUN(t,t+dt);
    file << t << "\t" << Rd.GetExpectationValue(Rd.d0,0,0) << "  " 
	 << Rd.GetExpectationValue(Rd.evol_b0_proj[0],0,0) << "  " 
	 << Rd.GetExpectationValue(Rd.evol_b0_proj[1],0,0) << endl;
    Rd.set_evol();
  }
  file.close();

  
  //Ask for runnint the gnuplot script
  string plt;
  cout << endl <<  "Done! " << endl <<  "  Do you want to run the gnuplot script? yes/no" << endl;
  cin >> plt;
  if(plt=="yes" || plt=="y"){
    return system("./plot.plt");
  }else{
    return 0;
  }

}
