#include <iostream>
#include <SQuIDs/SUNalg.h>

int main(){
  int dim=6;
  int Ngenerators=dim*dim;
  for(unsigned int i=0;i<Ngenerators;i++){
    for(unsigned int j=0;j<Ngenerators;j++){
      for(unsigned int k=0;k<Ngenerators;k++){
	SU_vector v1=SU_vector::Generator(dim,i);
	SU_vector v2=SU_vector::Generator(dim,j);
	SU_vector v3=SU_vector::Generator(dim,k);
	std::cout << i << "\t" << j << "\t"<< k << "\t" << v1*iCommutator(v2,v3) << std::endl;
      }
    }
  }
	
  return(0);
}
