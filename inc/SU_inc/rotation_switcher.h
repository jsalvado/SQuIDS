switch (dim){
case 2:
switch (i){
case 1:
switch (j){
case 2:
#include "RotationSU2_12.txt"
break;
}
break;
}
break;
case 3:
switch (i){
case 1:
switch (j){
case 2:
#include "RotationSU3_12.txt"
break;
case 3:
#include "RotationSU3_13.txt"
break;
}
break;
case 2:
switch (j){
case 3:
#include "RotationSU3_23.txt"
break;
}
break;
}
break;
case 4:
switch (i){
case 1:
switch (j){
case 2:
#include "RotationSU4_12.txt"
break;
case 3:
#include "RotationSU4_13.txt"
break;
case 4:
#include "RotationSU4_14.txt"
break;
}
break;
case 2:
switch (j){
case 3:
#include "RotationSU4_23.txt"
break;
case 4:
#include "RotationSU4_24.txt"
break;
}
break;
case 3:
switch (j){
case 4:
#include "RotationSU4_34.txt"
break;
}
break;
}
break;
case 5:
switch (i){
case 1:
switch (j){
case 2:
#include "RotationSU5_12.txt"
break;
case 3:
#include "RotationSU5_13.txt"
break;
case 4:
#include "RotationSU5_14.txt"
break;
case 5:
#include "RotationSU5_15.txt"
break;
}
break;
case 2:
switch (j){
case 3:
#include "RotationSU5_23.txt"
break;
case 4:
#include "RotationSU5_24.txt"
break;
case 5:
#include "RotationSU5_25.txt"
break;
}
break;
case 3:
switch (j){
case 4:
#include "RotationSU5_34.txt"
break;
case 5:
#include "RotationSU5_35.txt"
break;
}
break;
case 4:
switch (j){
case 5:
#include "RotationSU5_45.txt"
break;
}
break;
}
break;
case 6:
switch (i){
case 1:
switch (j){
case 2:
#include "RotationSU6_12.txt"
break;
case 3:
#include "RotationSU6_13.txt"
break;
case 4:
#include "RotationSU6_14.txt"
break;
case 5:
#include "RotationSU6_15.txt"
break;
case 6:
#include "RotationSU6_16.txt"
break;
}
break;
case 2:
switch (j){
case 3:
#include "RotationSU6_23.txt"
break;
case 4:
#include "RotationSU6_24.txt"
break;
case 5:
#include "RotationSU6_25.txt"
break;
case 6:
#include "RotationSU6_26.txt"
break;
}
break;
case 3:
switch (j){
case 4:
#include "RotationSU6_34.txt"
break;
case 5:
#include "RotationSU6_35.txt"
break;
case 6:
#include "RotationSU6_36.txt"
break;
}
break;
case 4:
switch (j){
case 5:
#include "RotationSU6_45.txt"
break;
case 6:
#include "RotationSU6_46.txt"
break;
}
break;
case 5:
switch (j){
case 6:
#include "RotationSU6_56.txt"
break;
}
break;
}
break;
default: 
throw std::runtime_error("SUN_rotation error. \n");
break;
}
