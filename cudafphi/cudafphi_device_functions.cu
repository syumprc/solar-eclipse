

#include <cuda_runtime.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdio.h>
static void print_device_info(int device_id){
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device_id);
  printf(" Device ID: %d\n", device_id);
  printf(" Device name: %s\n", prop.name);
  printf(" Device architecture: %d.%d\n", prop.major, prop.minor);
  printf(" Device global memory: %d\n", prop.totalGlobalMem);
  printf(" Device clockrate: %d\n", prop.clockRate);
}
static std::vector<int> print_devices(){
  int devices;
  cudaGetDeviceCount(&devices);
  printf("CUDA Capable Devices with architecture greater than equal to 3.5\n\n");
  //int usable_devices = 0;
  std::vector<int> usable_devices;
  cudaDeviceProp prop;

  for(int device = 0; device < devices; device++){
      cudaGetDeviceProperties(&prop, device);

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
      for(int iter = 0 ; iter < usable_devices.size(); iter++){

	  print_device_info(usable_devices[iter]);
	  get_input = true;

	  do{
	      try{
		  printf("y/n? ");
		  std::cin >> device_selected;
		  printf("\n");
		  std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
		  if(device_selected == "Y" || device_selected == "YES"){
		      selected_devices.push_back(usable_devices[iter]);
		  }else if (device_selected != "N" && device_selected != "NO"){
		      throw;
		  }
		  get_input = false;
	      }catch(...){
		  printf("Error in selection\n");
	      }
	  }while(get_input);

      }

      get_input = true;
      printf("Devices selected:\n");
      do{
	  for(int iter = 0; iter < selected_devices.size(); iter++){
	      print_device_info(selected_devices[iter]);
	  }
	  try{
	      printf("Is this right (y/n)? ");
	      std::cin >> device_selected;
	      printf("\n");
	      std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
	      if(device_selected == "Y" || device_selected == "YES"){
		 devices_selected = true;
	      }else if (device_selected != "N" && device_selected != "NO"){
		  throw;
	      }
	      get_input = false;
	  }catch(...){
	      printf("Error in selection\n");
	  }
      }while(get_input);
      if(devices_selected == false)
	selected_devices.clear();
  }while(devices_selected == false);

  return selected_devices;
}
