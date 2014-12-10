#include <iostream>
#include <fstream>
#include <SQuIDs/SUNalg.h>


int fact2(int i, int N){
  if(i==N-3)
    return 1;
  return fact2(i-1, N)*i;
}

int main(){
  int dim;
  int Ngenerators=dim*dim;
  std::ifstream file("FNumber.txt");
  int i,j,k;
  double val;

  for(int idim=2; idim<7; idim++){
    file >> dim;
    Ngenerators=dim*dim;
    int per=fact2(Ngenerators-1,Ngenerators-1);
    std::cout <<"tal "<< per << "  " << dim << std::endl;
    if(dim!=idim)
      std::cout << "Problem with dimensions  " << idim << "  " << dim << std::endl;
    for(unsigned int fi=0; fi<per; fi++){ 
      file >> i >> j >> k;
      file >> val;
      SU_vector v1=SU_vector::Generator(dim,i);
      SU_vector v2=SU_vector::Generator(dim,j);
      SU_vector v3=SU_vector::Generator(dim,k);
      if(v1*iCommutator(v2,v3)-val > 1e-4)
	std::cout << i << "\t" << j << "\t"<< k << "\t" << iCommutator(v2,v3)*v1-val  << std::endl;
    } 
  }
	
  return(0);
}
