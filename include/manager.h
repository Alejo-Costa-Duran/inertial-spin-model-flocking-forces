#pragma once

#include <SFML/Graphics.hpp>
#include <string>
#include <map>

class Manager {
public:
    Manager(const std::string& configFile = "config.txt", int windowWidth = 800, int windowHeight = 800);
    void Run();

private:
    int windowWidth_;
    int windowHeight_;
    double boxSize_;
    std::map<std::string, double> parameters_;

    void readConfig(const std::string& file);
    sf::Vector2f mapToScreen(double x, double y) const;
};