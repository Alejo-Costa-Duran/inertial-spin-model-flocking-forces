#include "manager.h"
#include "flock.h"
#include "bird.h"
#include "forcefield.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

Manager::Manager(const std::string& configFile, int windowWidth, int windowHeight)
    : windowWidth_(windowWidth), windowHeight_(windowHeight), boxSize_(50.0) {
    readConfig(configFile);
    if (parameters_.count("BoxSize")) {
        boxSize_ = parameters_["BoxSize"];
    }
}

void Manager::readConfig(const std::string& file) {
    std::ifstream in(file);
    if (!in.is_open()) {
        std::cout << "Manager Warning: Could not open " << file << ". Using default parameters.\n";
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '/') continue;
        std::stringstream ss(line);
        std::string key;
        double value;
        if (ss >> key >> value) {
            parameters_[key] = value;
        }
    }
    in.close();
}

sf::Vector2f Manager::mapToScreen(double x, double y) const {
    float screenX = static_cast<float>(windowWidth_ * (x + boxSize_ / 2.0) / boxSize_);
    float screenY = static_cast<float>(windowHeight_ * (boxSize_ / 2.0 - y) / boxSize_);
    return {screenX, screenY};
}

void Manager::Run() {
    // SFML Window Setup
    sf::RenderWindow window(sf::VideoMode(windowWidth_, windowHeight_), "ISM Swirl Dynamics - Interactive Visualizer");
    window.setFramerateLimit(120);

    // Extract parameters from config map using parameters_
    unsigned long int seed = static_cast<unsigned long int>(parameters_.count("Seed") ? parameters_["Seed"] : 1337);
    int numBirds = static_cast<int>(parameters_.count("NumBirds") ? parameters_["NumBirds"] : 50);
    double dt = parameters_.count("TimeStep") ? parameters_["TimeStep"] : 0.001;
    double maxSpeed = parameters_.count("MaxSpeed") ? parameters_["MaxSpeed"] : 1.0;
    
    double inertia = parameters_.count("VisualInertia") ? parameters_["VisualInertia"] : 
                     (parameters_.count("InertiaInit") ? parameters_["InertiaInit"] : 1.0);
                     
    double temperature = parameters_.count("VisualTemp") ? parameters_["VisualTemp"] : 
                         (parameters_.count("TempInit") ? parameters_["TempInit"] : 0.1);

    double friction = parameters_.count("Friction") ? parameters_["Friction"] : 1.0;
    double coupling = parameters_.count("Interaction") ? parameters_["Interaction"] : 1.0;

    // Generalized Well Parameters from config
    double wellStrength = parameters_.count("WellStrength") ? parameters_["WellStrength"] : 100.0;
    double wellRadius = parameters_.count("WellRadius") ? parameters_["WellRadius"] : 10.0;
    int wellExponent = parameters_.count("WellExponent") ? static_cast<int>(parameters_["WellExponent"]) : 12;

    std::cout << "Visualizer running with Generalized Well -> Strength: " << wellStrength 
              << ", Radius (R_0): " << wellRadius 
              << ", Exponent (n): " << wellExponent << "\n";

    // Initialize Flock
    Flock flock(seed, boxSize_, dt, maxSpeed, inertia, friction, temperature, coupling, 1.0, 0.0);

    // Spawn birds safely inside the well radius
    for (int idx = 0; idx < numBirds; ++idx) {
        double r0 = 2.0 + (gsl_rng_uniform(flock.m_mt) * (wellRadius * 0.7)); 
        double theta0 = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
        Pvector pos(r0 * std::cos(theta0), r0 * std::sin(theta0));
        Pvector vel(-maxSpeed * std::sin(theta0), maxSpeed * std::cos(theta0));
        Pvector acc(-maxSpeed * maxSpeed * std::cos(theta0) / r0, -maxSpeed * maxSpeed * std::sin(theta0) / r0);
        
        flock.bandada.push_back(Bird(false, idx, pos, vel, acc, inertia, maxSpeed, boxSize_));
    }

    // Use Generalized Power Confinement
    auto compositeForce = std::make_shared<CompositeForce>();
    compositeForce->addForce(std::make_shared<GeneralizedPowerConfinement>(wellStrength, wellRadius, wellExponent));
    const ForceField* activeField = compositeForce.get();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) window.close();
        }

        // Advance physics steps per frame
        for (int i = 0; i < 100; ++i) {
            flock.updateFlock(dt, activeField);
        }

        // Render Graphics
        window.clear(sf::Color(30, 30, 30));

        // Draw boundary reference circle matching WellRadius (R_0)
        sf::Vector2f center = mapToScreen(0.0, 0.0);
        float radiusPix = static_cast<float>((wellRadius / boxSize_) * windowWidth_);
        
        sf::CircleShape boundary(radiusPix);
        boundary.setOrigin(radiusPix, radiusPix);
        boundary.setPosition(center);
        boundary.setFillColor(sf::Color::Transparent);
        boundary.setOutlineColor(sf::Color(0, 150, 255, 150)); 
        boundary.setOutlineThickness(2.0f);
        window.draw(boundary);

        // Render birds
        for (const auto& bird : flock.bandada) {
            sf::Vector2f pos = mapToScreen(bird.position.x, bird.position.y);

            sf::CircleShape birdShape(4.0f);
            birdShape.setOrigin(4.0f, 4.0f);
            birdShape.setPosition(pos);

            double spin = bird.getSpin(inertia, maxSpeed);
            if (spin >= 0.0) {
                birdShape.setFillColor(sf::Color(0, 188, 212)); // Cyan
            } else {
                birdShape.setFillColor(sf::Color(255, 87, 34));  // Orange
            }
            window.draw(birdShape);

            float headLen = 10.0f;
            double velAngle = std::atan2(bird.velocity.y, bird.velocity.x);
            sf::Vertex line[] = {
                sf::Vertex(pos, sf::Color::White),
                sf::Vertex(sf::Vector2f(
                    pos.x + headLen * static_cast<float>(std::cos(velAngle)),
                    pos.y - headLen * static_cast<float>(std::sin(velAngle))
                ), sf::Color(200, 200, 200, 150))
            };
            window.draw(line, 2, sf::Lines);
        }

        window.display();
    }
}