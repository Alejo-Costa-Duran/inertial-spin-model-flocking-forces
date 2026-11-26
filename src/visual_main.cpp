#include "manager.h"
#include <iostream>
#include <exception>

int main(int argc, char* argv[]) {
    std::string configFile = (argc > 1) ? argv[1] : "config.txt";
    std::cout << "Loading visualizer configuration from: " << configFile << "\n";
    
    try {
        Manager managerWindow(configFile);
        managerWindow.Run();
    } catch (const std::exception& error) {
        std::cerr << "Visualizer failed: " << error.what() << "\n";
        return 1;
    }
    return 0;
}