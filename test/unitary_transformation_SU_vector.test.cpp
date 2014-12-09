#include <iostream>
#include <SQuIDs/SUNalg.h>

int main(){
  int dim=6;
  int Ngenerators=dim*dim;
  double th=1.456;
  double dl=1.231;
  double t=3.456;

  for(unsigned int i=0;i<Ngenerators;i++){
    for(unsigned int j=0;j<Ngenerators;j++){
      SU_vector v1=SU_vector::Generator(dim,i);
      SU_vector v2=SU_vector::Generator(dim,j);
      double sprod=v1*v2;
      for(unsigned int it=0;it<dim;it++){
	SU_vector v2=SU_vector::Projector(dim,it);
	double out=sprod-v1.SUEvolve(v2,t)*v2.SUEvolve(v1,t);
	if(fabs(out)>1e-15){
	    std::cout << "Generators: " << i << "  " << j << "  Evolution Operator: " << it << 
	      "  \t" << out  << std::endl;
	  }else{
	    std::cout  << "Generators: " << i << "  " << j << "  Evolution Operator: " << it << 
	      "  \t OK" <<  std::endl;
	  }
      }
    }
  }


  for(unsigned int i=0;i<Ngenerators;i++){
    for(unsigned int j=0;j<Ngenerators;j++){
      SU_vector v1=SU_vector::Generator(dim,i);
      SU_vector v2=SU_vector::Generator(dim,j);
      double sprod=v1*v2;
      for(unsigned int ir=1;ir<Ngenerators;ir++){
	for(unsigned int jr=ir+1;jr<Ngenerators;jr++){
	  double out=sprod-v1.Rotate(ir,jr,th,dl)*v2.Rotate(ir,jr,th,dl);
	  if(fabs(out)>1e-15){
	    std::cout << "Generators: " << i << "  " << j << "  Rotation: " << ir << "  " << jr << 
	      "  \t" << out  << std::endl;
	  }else{
	    std::cout  << "Generators: " << i << "  " << j << "  Rotation: " << ir << "  " << jr << 
	      "  \t OK" << std::endl;
	  }
	}
      }
    }
  }
  
  return(0);
}
