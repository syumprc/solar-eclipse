

#include <cuda_runtime.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <ncurses.h>
#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#endif
static inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{

   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      fprintf(stdout,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      endwin();
      if (abort) exit(code);
   }else{
	   fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
   }
}

static void print_device_info(int device_id){
  cudaDeviceProp prop;
  gpuErrchk(cudaGetDeviceProperties(&prop, device_id));
  printw(" Device ID: %d\n", device_id);
  printw(" Device name: %s\n", prop.name);
  printw(" Device architecture: %u.%u\n", prop.major, prop.minor);
  printw(" Device global memory: %u\n", prop.totalGlobalMem);
  printw(" Device clockrate: %u\n", prop.clockRate);

}
static std::vector<int> print_devices(){
  int devices;
  gpuErrchk(cudaGetDeviceCount(&devices));
  printw("CUDA Capable Devices with architecture greater than equal to 3.5\n\n");
  std::vector<int> usable_devices;
  cudaDeviceProp prop;

  for(int device = 0; device < devices; device++){
      gpuErrchk(cudaGetDeviceProperties(&prop, device));

      if(prop.major > 3){
    	  usable_devices.push_back(device);
	  	  printw(" Device ID: %d\n", device);
	  	  printw(" Device name: %s\n", prop.name);
	  	  printw(" Device architecture: %i.%i\n", prop.major, prop.minor);
      }else if((prop.major == 3) && (prop.minor >= 5)){
    	  printw(" Device ID: %i\n", device);
    	  printw(" Device name: %i\n", prop.name);
          printw(" Device architecture: %i.%i\n", prop.major, prop.minor);
          usable_devices.push_back(device);
      }
  	printw("Hit enter to continue\n");

    getch();

  }
  printw("\n");
  return usable_devices;
}

extern "C" std::vector<int> select_devices(){
  std::vector<int> selected_devices;

  if(initscr() == NULL){
	  printw("Failed to initialize select screen\n");
	  return selected_devices;
  }
  raw();
  keypad(stdscr, TRUE);
  noecho();
  bool devices_selected = false;
  bool get_input = true;
  std::vector<int> usable_devices = print_devices();
  if(usable_devices.size() == 0)
    return selected_devices;
  std::string device_selected;
  printw("\n");
  do{
      printw("Select devices:\n");
      for(std::vector<int>::iterator iter = usable_devices.begin() ; iter != usable_devices.end(); iter++){
    	  print_device_info(*iter);
    	  get_input = true;
    	  do{
    		  try{
    			  printw("y/n?\n");
    			  char selection = getch();
    			  device_selected = selection;

    			  std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
    			  if(device_selected == "Y" || device_selected == "YES"){
    				  selected_devices.push_back(*iter);
    				  get_input = false;
    			  }else if (device_selected == "N" || device_selected == "NO"){
    				  get_input = false;
    			  }else if (device_selected == "E"){
    				  selected_devices.clear();
    				  endwin();
    				  return selected_devices;
    			  }else{
    				  throw 0;
    			  }

    		  }catch(...){
    			  printw("Error in selection\n");
    			  continue;
    		  }

    	      refresh();

    	  }while(get_input);

      }

      get_input = true;
      printw("Devices selected:\n");
      do{
    	  for(std::vector<int>::iterator iter = selected_devices.begin(); iter != selected_devices.end(); iter++){
    		  print_device_info(*iter);
    	  }
	  	  try{
	  		  printw("Is this right (y/n)? \n");
			  char selection = getch();
			  device_selected = selection;
	  		  printw("\n");
	  		  std::transform(device_selected.begin(), device_selected.end(),device_selected.begin(), ::toupper);
	  		  if(device_selected == "Y" || device_selected == "YES"){
	  			  devices_selected = true;
	  			  get_input = false;
	  		  }else if (device_selected == "N" || device_selected == "NO"){
	  			  get_input = false;
	  		  }else if (device_selected == "E"){
				  selected_devices.clear();
				  endwin();
				  return selected_devices;
			  }else{
	  			  throw 0;
	  		  }
	  	  }catch(...){
	  		  printw("Error in selection\n");
	  		  continue;
	  	  }
	  	  refresh();
      }while(get_input);
      if(devices_selected == false) selected_devices.clear();
  }while(devices_selected == false);
  printw("Hit enter to continue\n");

  getch();
  endwin();
  return selected_devices;
}

extern "C" std::vector<int> select_all_devices(){
	std::vector<int> all_devices;


	raw();
	keypad(stdscr, TRUE);
	noecho();

	if(initscr() == NULL){
		printw("Failed to initialize select screen\n");

		return all_devices;
	}
	all_devices = print_devices();
	printw("Hit enter to continue\n");
	getch();
	endwin();
	return all_devices;
}
