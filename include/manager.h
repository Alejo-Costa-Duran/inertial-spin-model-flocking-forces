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
    double viewRadius_;
    std::map<std::string, double> parameters_;
    static constexpr int parameterPanelWidth_ = 300;

    void readConfig(const std::string& file);
    double pixelsPerUnit() const;
    sf::Vector2f mapToScreen(double x, double y) const;
};