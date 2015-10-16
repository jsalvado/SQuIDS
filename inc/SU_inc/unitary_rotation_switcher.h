switch (dim){
case 2:
#include "SUUnitaryRotate2.txt"
break;
case 3:
#include "SUUnitaryRotate3.txt"
break;
case 4:
#include "SUUnitaryRotate4.txt"
break;
case 5:
#include "SUUnitaryRotate5.txt"
break;
case 6:
#include "SUUnitaryRotate6.txt"
break;
default: 
throw std::runtime_error("SUN_unitary_rotation error. \n");
break;
}
