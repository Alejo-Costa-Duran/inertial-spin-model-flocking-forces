#include "manager.h"
#include <iostream>

int main(int argc, char* argv[]) {
    std::string configFile = (argc > 1) ? argv[1] : "config.txt";
    std::cout << "Loading visualizer configuration from: " << configFile << "\n";
    
    Manager managerWindow(configFile);
    managerWindow.Run();
    return 0;
}