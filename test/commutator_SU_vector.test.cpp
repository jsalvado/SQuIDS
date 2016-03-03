#include <iostream>
#include <fstream>
#include <SQuIDS/SUNalg.h>

int main(){
  using squids::SU_vector;
  
  int dim;
  int Ngenerators;
  std::ifstream file("FNumber.txt");
  int in,jn,kn;
  double val;
  bool ok;
  bool set=false;

  for(int d=2; d<7;d++){
    file >> dim; 
    //    std::cout << dim << std::endl;
    Ngenerators=dim*dim;
    ok=true;
    for(int i=0; i<Ngenerators; i++){
      SU_vector v1=SU_vector::Generator(dim,i);
      for(int j=0; j<Ngenerators; j++){
	SU_vector v2=SU_vector::Generator(dim,j);
	for(int k=0; k<Ngenerators; k++){
	  SU_vector v3=SU_vector::Generator(dim,k); 
	  //	std::cout << i << "  " << j << "  " << k << "  " << (iCommutator(v2,v3)*v1)/4.0 << std::endl;
	  if(std::abs(iCommutator(v1,v2)*v3)>1e-15){
	    if(!set){
	      file >> in; 
	      file >> jn; 
	      file >> kn;	    
	      file >> val;
	    }
	    if( i!=in || k!=kn ||k!=kn){
	      ok=false;
	    }else{
	      if(val+(iCommutator(v1,v2)*v3)/4.0>1e-4){
		ok=false;
		std::cout << i << "  " << j << "  " << k << "  " << val+(iCommutator(v2,v3)*v1)/4.0 << "  " << val << std::endl;
	      }else{
		set=false;
	      }
	    }
	  }

	}
      }
    }
  if(!ok)
    std::cout << "commutator wrong, dim: " << dim << std::endl;
  
  }

  return(0);
}
