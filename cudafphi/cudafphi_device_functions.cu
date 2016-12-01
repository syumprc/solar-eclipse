

#include <cuda_runtime.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "cudafphi.cuh"
static void print_device_info(int device_id){
  cudaDeviceProp prop;
  gpuErrchk(cudaGetDeviceProperties(&prop, device_id));
  printf(" Device ID: %d\n", device_id);
  printf(" Device name: %s\n", prop.name);
  printf(" Device architecture: %d.%d\n", prop.major, prop.minor);
  printf(" Device global memory: %d\n", prop.totalGlobalMem);
  printf(" Device clockrate: %d\n", prop.clockRate);
}
static std::vector<int> print_devices(){
  int devices;
  gpuErrchk(cudaGetDeviceCount(&devices));
  printf("CUDA Capable Devices with architecture greater than equal to 3.5\n\n");
  std::vector<int> usable_devices;
  cudaDeviceProp prop;

  for(int device = 0; device < devices; device++){
      gpuErrchk(cudaGetDeviceProperties(&prop, device));

      if(prop.major > 3){
	  usable_devices.push_back(device);
	  printf(" Device ID: %d\n", device);
	  printf(" Device name: %s\n", prop.name);
	  printf(" Device architecture: %d.%d\n", prop.major, prop.minor);
      }else if((prop.major == 3) && (prop.minor >= 5)){
	  printf(" Device ID: %d\n", device);
	  printf(" Device name: %s\n", prop.name);
          printf(" Device architecture: %d.%d\n", prop.major, prop.minor);
	  usable_devices.push_back(device);
      }

  }
  printf("\n");
  return usable_devices;
}

extern "C" std::vector<int> select_devices(){
  bool devices_selected = false;
  bool get_input = true;
  std::vector<int> usable_devices = print_devices();
  std::vector<int> selected_devices;
  if(usable_devices.size() == 0)
    return selected_devices;
  std::string device_selected;
  printf("\n");
  do{
      printf("Select devices:\n");
      for(std::vector<int>::iterator iter = usable_devices.begin() ; iter != usable_devices.end(); iter++){
    	  print_device_info(*iter);
    	  get_input = true;
    	  do{
    		  try{
    			  printf("y/n? ");
    			  std::cin >> device_selected;
    			  printf("\n");
    			  std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
    			  if(device_selected == "Y" || device_selected == "YES"){
    				  selected_devices.push_back(*iter);
    				  get_input = false;
    			  }else if (device_selected == "N" || device_selected == "NO"){
    				  get_input = false;
    			  }else{
    				  throw 0;
    			  }

    		  }catch(...){
    			  printf("Error in selection\n");
    			  continue;
    		  }
    	  }while(get_input);

      }

      get_input = true;
      printf("Devices selected:\n");
      do{
    	  for(std::vector<int>::iterator iter = selected_devices.begin(); iter != selected_devices.end(); iter++){
    		  print_device_info(*iter);
    	  }
	  	  try{
	  		  printf("Is this right (y/n)? ");
	  		  std::cin >> device_selected;
	  		  printf("\n");
	  		  std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
	  		  if(device_selected == "Y" || device_selected == "YES"){
	  			  devices_selected = true;
	  			  get_input = false;
	  		  }else if (device_selected == "N" || device_selected == "NO"){
	  			  get_input = false;
	  		  }else{
	  			  throw 0;
	  		  }
	  	  }catch(...){
	  		  printf("Error in selection\n");
	  		  continue;
	  	  }
      }while(get_input);
      if(devices_selected == false) selected_devices.clear();
  }while(devices_selected == false);

  return selected_devices;
}

extern "C" std::vector<int> select_all_devices(){
	std::cout << "Using the following devices:\n";
	return print_devices();
}
