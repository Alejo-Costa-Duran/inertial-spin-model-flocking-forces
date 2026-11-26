#include "manager.h"
#include "flock.h"
#include "bird.h"
#include "forcefield.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>

Manager::Manager(const std::string& configFile, int windowWidth, int windowHeight)
    : windowWidth_(windowWidth), windowHeight_(windowHeight), boxSize_(50.0), viewRadius_(50.0) {
    readConfig(configFile);
    if (parameters_.count("BoxSize")) {
        boxSize_ = parameters_["BoxSize"];
    }
    viewRadius_ = boxSize_;
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

double Manager::pixelsPerUnit() const {
    constexpr double padding = 12.0;
    return std::min(windowWidth_ - 2.0 * padding, windowHeight_ - 2.0 * padding) /
           (2.0 * viewRadius_);
}

sf::Vector2f Manager::mapToScreen(double x, double y) const {
    const double scale = pixelsPerUnit();
    const double boxPixels = 2.0 * viewRadius_ * scale;
    const double offsetX = (windowWidth_ - boxPixels) / 2.0;
    const double offsetY = (windowHeight_ - boxPixels) / 2.0;
    float screenX = static_cast<float>(offsetX + (x + viewRadius_) * scale);
    float screenY = static_cast<float>(offsetY + (viewRadius_ - y) * scale);
    return {screenX, screenY};
}

void Manager::Run() {
    constexpr int stepsPerFrame = 100;

    sf::Font font;
    const std::vector<std::string> fontPaths = {
        "C:\\Windows\\Fonts\\arial.ttf",
        "C:\\Windows\\Fonts\\segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf"
    };
    bool fontLoaded = false;
    for (const auto& fontPath : fontPaths) {
        if (font.loadFromFile(fontPath)) {
            fontLoaded = true;
            break;
        }
    }
    if (!fontLoaded) {
        throw std::runtime_error("Visualizer could not load a system font for the parameter panel.");
    }

    // SFML Window Setup
    sf::RenderWindow window(
        sf::VideoMode(windowWidth_ + parameterPanelWidth_, windowHeight_),
        "ISM Swirl Dynamics - Interactive Visualizer"
    );
    window.setFramerateLimit(120);

    // The visualizer uses dimensionless units with v0 = J = eta = 1.
    unsigned long int seed = static_cast<unsigned long int>(parameters_.count("Seed") ? parameters_["Seed"] : 1337);
    int numBirds = static_cast<int>(parameters_.count("NumBirds") ? parameters_["NumBirds"] : 50);
    double dt = parameters_.count("TimeStep") ? parameters_["TimeStep"] : 0.001;
    constexpr double maxSpeed = 1.0;
    
    double omega = parameters_.count("VisualOmega") ? parameters_["VisualOmega"] :
                   (parameters_.count("Omega") ? parameters_["Omega"] :
                   (parameters_.count("VisualInertia") ? parameters_["VisualInertia"] :
                   (parameters_.count("InertiaInit") ? parameters_["InertiaInit"] : 1.0)));
                     
    double temperature = parameters_.count("VisualTemp") ? parameters_["VisualTemp"] : 
                         (parameters_.count("TempInit") ? parameters_["TempInit"] : 0.1);

    constexpr double friction = 1.0;
    constexpr double coupling = 1.0;
    const bool periodicBoundary = parameters_.count("PeriodicBoundary") &&
                                  parameters_["PeriodicBoundary"] != 0.0;
    const bool randomInitialConditions = parameters_.count("RandomInitialConditions") &&
                                         parameters_["RandomInitialConditions"] != 0.0;

    // Generalized Well Parameters from config
    double wellStrength = parameters_.count("WellStrength") ? parameters_["WellStrength"] : 100.0;
    double wellRadius = parameters_.count("WellRadius") ? parameters_["WellRadius"] : 10.0;
    double wellExponent = static_cast<int>(parameters_.count("WellExponent") ? parameters_["WellExponent"] : 12);

    if (periodicBoundary && wellStrength > 0.0) {
        throw std::invalid_argument(
            "PeriodicBoundary=1 cannot be combined with an active generalized well; set WellStrength=0."
        );
    }

    std::cout << "Visualizer units: v0 = J = eta = 1; Omega: " << omega
              << ", T/J: " << temperature << "\n";
    std::cout << "Boundary mode: " << (periodicBoundary ? "periodic" : "open") << "\n";
    std::cout << "Generalized Well -> Strength*: " << wellStrength
              << ", Radius r0*: " << wellRadius
              << ", Exponent (n): " << static_cast<int>(wellExponent) << "\n";

    // Initialize Flock
    Flock flock(seed, boxSize_, dt, maxSpeed, omega, friction, temperature, coupling, 1.0, 0.0, periodicBoundary);

    auto makeForceField = [&]() {
        auto forceField = std::make_shared<CompositeForce>();
        if (wellStrength > 0.0) {
            forceField->addForce(std::make_shared<GeneralizedPowerConfinement>(
                wellStrength, wellRadius, static_cast<int>(wellExponent)
            ));
        }
        return forceField;
    };
    std::shared_ptr<CompositeForce> compositeForce = makeForceField();
    bool showNeighborRadius = false;
    int inspectedBirdIndex = -1;

    auto respawnFlock = [&]() {
        flock.bandada.clear();
        gsl_rng_set(flock.m_mt, seed);
        for (int idx = 0; idx < numBirds; ++idx) {
            Pvector pos(0.0, 0.0);
            Pvector vel(0.0, 0.0);
            Pvector acc(0.0, 0.0);
            if (randomInitialConditions) {
                const double initialRadius = std::min(wellRadius, boxSize_);
                const double radialSample = initialRadius * std::sqrt(gsl_rng_uniform(flock.m_mt));
                const double positionAngle = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
                pos.x = radialSample * std::cos(positionAngle);
                pos.y = radialSample * std::sin(positionAngle);
                const double heading = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
                vel.x = maxSpeed * std::cos(heading);
                vel.y = maxSpeed * std::sin(heading);
            } else {
                const double r0 = 2.0 + (gsl_rng_uniform(flock.m_mt) * (wellRadius * 0.7));
                const double theta0 = 2.0 * M_PI * gsl_rng_uniform(flock.m_mt);
                pos.x = r0 * std::cos(theta0);
                pos.y = r0 * std::sin(theta0);
                vel.x = -maxSpeed * std::sin(theta0);
                vel.y = maxSpeed * std::cos(theta0);
                acc.x = -maxSpeed * maxSpeed * std::cos(theta0) / r0;
                acc.y = -maxSpeed * maxSpeed * std::sin(theta0) / r0;
            }
            flock.bandada.push_back(Bird(false, idx, pos, vel, acc));
        }
        if (periodicBoundary) flock.boundary();
        flock.calculateForces(compositeForce.get());
        if (randomInitialConditions) flock.initializeAccelerations();
        inspectedBirdIndex = flock.bandada.empty() ? -1 : 0;
    };
    respawnFlock();

    struct SliderControl {
        std::string label;
        double* value;
        double minimum;
        double maximum;
        bool logarithmic;
        bool integer;
        bool updatesFlockConstants;
        bool updatesForceField;
    };
    std::vector<SliderControl> sliders = {
        {"Inertia Omega (log)", &omega, 0.0001, 1.0, true, false, true, false},
        {"Temperature T/J", &temperature, 0.0, 3.0, false, false, true, false}
    };
    if (!periodicBoundary) {
        sliders.push_back({"Well radius r0*", &wellRadius, 0.1, 15.0, false, false, false, true});
        sliders.push_back({"Well exponent", &wellExponent, 1.0, 20.0, false, true, false, true});
        sliders.push_back({"Well strength s*", &wellStrength, 0.0, 200.0, false, false, false, true});
    }

    const float sliderLeft = static_cast<float>(windowWidth_) + 20.0f;
    const float sliderWidth = static_cast<float>(parameterPanelWidth_) - 40.0f;
    const float sliderRowHeight = 68.0f;
    const float sliderTop = 164.0f;
    const float resetButtonX = static_cast<float>(windowWidth_) + 205.0f;
    const float resetButtonY = 20.0f;
    const float neighborButtonX = static_cast<float>(windowWidth_) + 20.0f;
    const float neighborButtonY = 500.0f;
    const float neighborButtonWidth = 150.0f;
    const float neighborButtonHeight = 32.0f;
    int activeSlider = -1;

    auto formatValue = [](double value) {
        std::ostringstream formatted;
        formatted << std::setprecision(5) << value;
        return formatted.str();
    };

    auto setSliderFromMouse = [&](std::size_t index, int mouseX) {
        SliderControl& slider = sliders[index];
        double position = (static_cast<double>(mouseX) - sliderLeft) / sliderWidth;
        if (position < 0.0) position = 0.0;
        if (position > 1.0) position = 1.0;

        double value = slider.logarithmic
            ? slider.minimum * std::pow(slider.maximum / slider.minimum, position)
            : slider.minimum + position * (slider.maximum - slider.minimum);
        if (slider.integer) value = std::round(value);
        if (value == *slider.value) return;
        *slider.value = value;

        if (slider.updatesFlockConstants) {
            flock.inertia = omega;
            flock.temperature = temperature;
            flock.setConstants(dt);
        }
        if (slider.updatesForceField) compositeForce = makeForceField();
    };

    sf::Clock fpsClock;
    unsigned int fpsFrameCount = 0;
    float currentFps = 0.0f;
    std::string simulationError;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) window.close();
            if (event.type == sf::Event::MouseButtonReleased &&
                event.mouseButton.button == sf::Mouse::Left) {
                activeSlider = -1;
            }
            if (event.type == sf::Event::MouseButtonPressed &&
                event.mouseButton.button == sf::Mouse::Left) {
                activeSlider = -1;
                if (event.mouseButton.x >= resetButtonX &&
                    event.mouseButton.x <= resetButtonX + 75.0f &&
                    event.mouseButton.y >= resetButtonY &&
                    event.mouseButton.y <= resetButtonY + 32.0f) {
                    respawnFlock();
                    simulationError.clear();
                    activeSlider = -1;
                } else if (event.mouseButton.x >= neighborButtonX &&
                           event.mouseButton.x <= neighborButtonX + neighborButtonWidth &&
                           event.mouseButton.y >= neighborButtonY &&
                           event.mouseButton.y <= neighborButtonY + neighborButtonHeight) {
                    showNeighborRadius = !showNeighborRadius;
                } else {
                    bool sliderHit = false;
                    for (std::size_t index = 0; index < sliders.size(); ++index) {
                        const float rowY = sliderTop + static_cast<float>(index) * sliderRowHeight;
                        if (event.mouseButton.x >= sliderLeft &&
                            event.mouseButton.x <= sliderLeft + sliderWidth &&
                            event.mouseButton.y >= rowY + 20.0f &&
                            event.mouseButton.y <= rowY + 38.0f) {
                            activeSlider = static_cast<int>(index);
                            setSliderFromMouse(index, event.mouseButton.x);
                            sliderHit = true;
                            break;
                        }
                    }
                    if (!sliderHit && showNeighborRadius && event.mouseButton.x < windowWidth_) {
                        const float selectionRadius = 12.0f;
                        float nearestDistanceSquared = selectionRadius * selectionRadius;
                        for (std::size_t index = 0; index < flock.bandada.size(); ++index) {
                            const sf::Vector2f position = mapToScreen(
                                flock.bandada[index].position.x,
                                flock.bandada[index].position.y
                            );
                            const float dx = position.x - static_cast<float>(event.mouseButton.x);
                            const float dy = position.y - static_cast<float>(event.mouseButton.y);
                            const float distanceSquared = dx * dx + dy * dy;
                            if (distanceSquared <= nearestDistanceSquared) {
                                nearestDistanceSquared = distanceSquared;
                                inspectedBirdIndex = static_cast<int>(index);
                            }
                        }
                    }
                }
            }
            if (event.type == sf::Event::MouseMoved && activeSlider >= 0) {
                setSliderFromMouse(static_cast<std::size_t>(activeSlider), event.mouseMove.x);
            }
        }

        // Advance physics steps per frame, but keep the visualizer open if the
        // configured forces and time step cannot satisfy the speed constraint.
        if (simulationError.empty()) {
            for (int i = 0; i < stepsPerFrame; ++i) {
                try {
                    flock.updateFlock(dt, compositeForce.get());
                } catch (const std::exception& error) {
                    simulationError = error.what();
                    std::cerr << "Visualizer paused: " << simulationError << "\n";
                    break;
                }
            }
        }

        if (periodicBoundary) {
            viewRadius_ = boxSize_;
        } else {
            double requiredViewRadius = boxSize_;
            for (const Bird& bird : flock.bandada) {
                requiredViewRadius = std::max(
                    requiredViewRadius,
                    std::max(std::abs(bird.position.x), std::abs(bird.position.y)) + 1.0
                );
            }
            if (requiredViewRadius > viewRadius_) {
                viewRadius_ = requiredViewRadius;
            } else {
                viewRadius_ = std::max(
                    requiredViewRadius,
                    viewRadius_ * 0.96 + requiredViewRadius * 0.04
                );
            }
        }

        std::vector<Bird> inspectedNeighbors;
        if (showNeighborRadius &&
            inspectedBirdIndex >= 0 &&
            static_cast<std::size_t>(inspectedBirdIndex) < flock.bandada.size()) {
            inspectedNeighbors = flock.flocking(inspectedBirdIndex);
        }

        // Render Graphics
        window.clear(sf::Color(30, 30, 30));

        const float boxPixels = static_cast<float>(2.0 * boxSize_ * pixelsPerUnit());
        sf::RectangleShape boxBoundary(sf::Vector2f(boxPixels, boxPixels));
        boxBoundary.setPosition(mapToScreen(-boxSize_, boxSize_));
        boxBoundary.setFillColor(sf::Color::Transparent);
        boxBoundary.setOutlineColor(sf::Color(90, 100, 115));
        boxBoundary.setOutlineThickness(1.0f);
        window.draw(boxBoundary);

        // Draw boundary reference circle matching WellRadius (R_0)
        sf::Vector2f center = mapToScreen(0.0, 0.0);
        float radiusPix = static_cast<float>(wellRadius * pixelsPerUnit());
        
        sf::CircleShape boundary(radiusPix);
        boundary.setOrigin(radiusPix, radiusPix);
        boundary.setPosition(center);
        boundary.setFillColor(sf::Color::Transparent);
        boundary.setOutlineColor(sf::Color(0, 150, 255, 150)); 
        boundary.setOutlineThickness(2.0f);
        window.draw(boundary);

        if (showNeighborRadius &&
            inspectedBirdIndex >= 0 &&
            static_cast<std::size_t>(inspectedBirdIndex) < flock.bandada.size()) {
            const Bird& inspectedBird = flock.bandada[inspectedBirdIndex];
            const double radius = flock.numNeigh;
            const float radiusPixels = static_cast<float>(radius * pixelsPerUnit());
            sf::CircleShape neighborRadius(radiusPixels);
            neighborRadius.setOrigin(radiusPixels, radiusPixels);
            neighborRadius.setFillColor(sf::Color::Transparent);
            neighborRadius.setOutlineColor(sf::Color(210, 210, 210, 190));
            neighborRadius.setOutlineThickness(1.0f);

            const int imageRange = periodicBoundary ? 1 : 0;
            const double period = 2.0 * boxSize_;
            for (int imageX = -imageRange; imageX <= imageRange; ++imageX) {
                for (int imageY = -imageRange; imageY <= imageRange; ++imageY) {
                    const double centerX = inspectedBird.position.x + imageX * period;
                    const double centerY = inspectedBird.position.y + imageY * period;
                    if (periodicBoundary &&
                        (centerX + radius < -boxSize_ || centerX - radius > boxSize_ ||
                         centerY + radius < -boxSize_ || centerY - radius > boxSize_)) {
                        continue;
                    }
                    neighborRadius.setPosition(mapToScreen(centerX, centerY));
                    window.draw(neighborRadius);
                }
            }
        }

        // Render birds
        for (const auto& bird : flock.bandada) {
            sf::Vector2f pos = mapToScreen(bird.position.x, bird.position.y);

            sf::CircleShape birdShape(4.0f);
            birdShape.setOrigin(4.0f, 4.0f);
            birdShape.setPosition(pos);

            double spin = bird.getSpin(omega, maxSpeed);
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

            if (showNeighborRadius) {
                const bool isInspectedBird =
                    inspectedBirdIndex >= 0 &&
                    static_cast<std::size_t>(inspectedBirdIndex) < flock.bandada.size() &&
                    bird.idx == flock.bandada[inspectedBirdIndex].idx;
                const bool isNeighbor = std::any_of(
                    inspectedNeighbors.begin(),
                    inspectedNeighbors.end(),
                    [&](const Bird& neighbor) { return neighbor.idx == bird.idx; }
                );
                if (isInspectedBird || isNeighbor) {
                    sf::CircleShape marker(isInspectedBird ? 6.5f : 7.5f);
                    const float markerRadius = isInspectedBird ? 6.5f : 7.5f;
                    marker.setOrigin(markerRadius, markerRadius);
                    marker.setPosition(pos);
                    marker.setFillColor(sf::Color::Transparent);
                    marker.setOutlineColor(sf::Color(220, 220, 220, 230));
                    marker.setOutlineThickness(1.5f);
                    window.draw(marker);
                }
            }
        }

        // Keep parameter text in a dedicated sidebar, outside the simulation area.
        const float panelX = static_cast<float>(windowWidth_);
        sf::RectangleShape panel(sf::Vector2f(
            static_cast<float>(parameterPanelWidth_),
            static_cast<float>(windowHeight_)
        ));
        panel.setPosition(panelX, 0.0f);
        panel.setFillColor(sf::Color(22, 26, 34));
        window.draw(panel);

        sf::Text heading("Dimensionless", font, 20);
        heading.setPosition(panelX + 20.0f, 22.0f);
        heading.setFillColor(sf::Color(230, 235, 245));
        window.draw(heading);

        sf::RectangleShape resetButton(sf::Vector2f(75.0f, 32.0f));
        resetButton.setPosition(resetButtonX, resetButtonY);
        resetButton.setFillColor(sf::Color(45, 93, 125));
        window.draw(resetButton);

        sf::Text resetButtonLabel("Reset flock", font, 13);
        resetButtonLabel.setPosition(resetButtonX + 7.0f, resetButtonY + 8.0f);
        resetButtonLabel.setFillColor(sf::Color::White);
        window.draw(resetButtonLabel);

        sf::Text subheading("Fixed scales: v0 = J = eta = 1", font, 12);
        subheading.setPosition(panelX + 20.0f, 54.0f);
        subheading.setFillColor(sf::Color(150, 165, 185));
        window.draw(subheading);

        sf::Text controlScales("Omega = J chi/eta^2; T* = T/J", font, 11);
        controlScales.setPosition(panelX + 20.0f, 70.0f);
        controlScales.setFillColor(sf::Color(150, 165, 185));
        window.draw(controlScales);

        sf::Text wellScales("r0* = J r0/(v0 eta); s* = v0 s/eta", font, 10);
        wellScales.setPosition(panelX + 20.0f, 85.0f);
        wellScales.setFillColor(sf::Color(150, 165, 185));
        window.draw(wellScales);

        sf::Text unitScales("t0 = eta/J = 1; L0 = v0 eta/J = 1", font, 10);
        unitScales.setPosition(panelX + 20.0f, 99.0f);
        unitScales.setFillColor(sf::Color(150, 165, 185));
        window.draw(unitScales);

        sf::Text fpsLabel("FPS", font, 15);
        fpsLabel.setPosition(panelX + 20.0f, 116.0f);
        fpsLabel.setFillColor(sf::Color(180, 192, 210));
        window.draw(fpsLabel);

        sf::Text fpsValue(std::to_string(static_cast<int>(currentFps + 0.5f)), font, 15);
        fpsValue.setPosition(panelX + 160.0f, 116.0f);
        fpsValue.setFillColor(sf::Color(245, 245, 245));
        window.draw(fpsValue);

        sf::Text boundaryMode(
            periodicBoundary ? "Boundary: periodic (well disabled)" : "Boundary: open",
            font,
            11
        );
        boundaryMode.setPosition(panelX + 20.0f, 138.0f);
        boundaryMode.setFillColor(sf::Color(150, 165, 185));
        window.draw(boundaryMode);

        for (std::size_t index = 0; index < sliders.size(); ++index) {
            const SliderControl& slider = sliders[index];
            const float rowY = sliderTop + static_cast<float>(index) * sliderRowHeight;
            double position;
            if (slider.logarithmic) {
                const double boundedValue = std::max(slider.minimum, std::min(slider.maximum, *slider.value));
                position = std::log(boundedValue / slider.minimum) /
                           std::log(slider.maximum / slider.minimum);
            } else {
                position = (*slider.value - slider.minimum) / (slider.maximum - slider.minimum);
                position = std::max(0.0, std::min(1.0, position));
            }

            sf::Text label(slider.label, font, 14);
            label.setPosition(sliderLeft, rowY);
            label.setFillColor(sf::Color(180, 192, 210));
            window.draw(label);

            sf::Text displayedValue(formatValue(*slider.value), font, 14);
            displayedValue.setPosition(panelX + 180.0f, rowY);
            displayedValue.setFillColor(sf::Color(245, 245, 245));
            window.draw(displayedValue);

            sf::RectangleShape track(sf::Vector2f(sliderWidth, 4.0f));
            track.setPosition(sliderLeft, rowY + 27.0f);
            track.setFillColor(sf::Color(65, 72, 84));
            window.draw(track);

            sf::RectangleShape filledTrack(sf::Vector2f(
                sliderWidth * static_cast<float>(position), 4.0f
            ));
            filledTrack.setPosition(sliderLeft, rowY + 27.0f);
            filledTrack.setFillColor(sf::Color(0, 150, 200));
            window.draw(filledTrack);

            sf::CircleShape handle(7.0f);
            handle.setOrigin(7.0f, 7.0f);
            handle.setPosition(sliderLeft + sliderWidth * static_cast<float>(position), rowY + 29.0f);
            handle.setFillColor(activeSlider == static_cast<int>(index)
                ? sf::Color(255, 190, 70)
                : sf::Color(100, 220, 245));
            window.draw(handle);

            sf::Text minimum(formatValue(slider.minimum), font, 10);
            minimum.setPosition(sliderLeft, rowY + 36.0f);
            minimum.setFillColor(sf::Color(130, 145, 165));
            window.draw(minimum);

            sf::Text maximum(formatValue(slider.maximum), font, 10);
            maximum.setPosition(
                sliderLeft + sliderWidth - maximum.getLocalBounds().width,
                rowY + 36.0f
            );
            maximum.setFillColor(sf::Color(130, 145, 165));
            window.draw(maximum);
        }

        sf::RectangleShape neighborButton(
            sf::Vector2f(neighborButtonWidth, neighborButtonHeight)
        );
        neighborButton.setPosition(neighborButtonX, neighborButtonY);
        neighborButton.setFillColor(showNeighborRadius
            ? sf::Color(78, 105, 118)
            : sf::Color(45, 93, 125));
        window.draw(neighborButton);

        sf::Text neighborButtonLabel(
            showNeighborRadius ? "Radius: ON" : "Radius: OFF",
            font,
            13
        );
        neighborButtonLabel.setPosition(neighborButtonX + 12.0f, neighborButtonY + 8.0f);
        neighborButtonLabel.setFillColor(sf::Color::White);
        window.draw(neighborButtonLabel);

        sf::Text neighborHint("Click a bird to inspect its neighbors", font, 10);
        neighborHint.setPosition(panelX + 20.0f, neighborButtonY + 38.0f);
        neighborHint.setFillColor(sf::Color(150, 165, 185));
        window.draw(neighborHint);

        sf::Text neighborCutoff("Cutoff: 1.0 simulation unit", font, 10);
        neighborCutoff.setPosition(panelX + 20.0f, neighborButtonY + 52.0f);
        neighborCutoff.setFillColor(sf::Color(150, 165, 185));
        window.draw(neighborCutoff);

        if (showNeighborRadius &&
            inspectedBirdIndex >= 0 &&
            static_cast<std::size_t>(inspectedBirdIndex) < flock.bandada.size()) {
            sf::Text neighborCount(
                "Bird " + std::to_string(flock.bandada[inspectedBirdIndex].idx) +
                    " sees " + std::to_string(inspectedNeighbors.size()) + " of " +
                    std::to_string(flock.bandada.size() - 1),
                font,
                10
            );
            neighborCount.setPosition(panelX + 20.0f, neighborButtonY + 66.0f);
            neighborCount.setFillColor(sf::Color(210, 210, 210));
            window.draw(neighborCount);
        }

        if (!simulationError.empty()) {
            sf::Text pausedMessage(
                "PAUSED: adjust settings, then Reset flock",
                font,
                10
            );
            pausedMessage.setPosition(panelX + 20.0f, 580.0f);
            pausedMessage.setFillColor(sf::Color(255, 150, 100));
            window.draw(pausedMessage);
        }

        sf::Text otherValuesHeading("Other active values", font, 14);
        otherValuesHeading.setPosition(panelX + 20.0f, 606.0f);
        otherValuesHeading.setFillColor(sf::Color(180, 192, 210));
        window.draw(otherValuesHeading);

        const std::vector<std::pair<std::string, double>> otherValues = {
            {"Seed", static_cast<double>(seed)},
            {"Birds", static_cast<double>(numBirds)},
            {"Box B*", boxSize_},
            {"Time step*", dt},
            {"v0=J=eta", 1.0},
            {"Steps/frame", static_cast<double>(stepsPerFrame)}
        };
        for (std::size_t index = 0; index < otherValues.size(); index += 2) {
            const float rowY = 632.0f + static_cast<float>(index / 2) * 28.0f;

            sf::Text firstLabel(otherValues[index].first, font, 11);
            firstLabel.setPosition(panelX + 20.0f, rowY);
            firstLabel.setFillColor(sf::Color(150, 165, 185));
            window.draw(firstLabel);

            sf::Text firstValue(formatValue(otherValues[index].second), font, 11);
            firstValue.setPosition(panelX + 82.0f, rowY);
            firstValue.setFillColor(sf::Color(220, 225, 235));
            window.draw(firstValue);

            sf::Text secondLabel(otherValues[index + 1].first, font, 11);
            secondLabel.setPosition(panelX + 155.0f, rowY);
            secondLabel.setFillColor(sf::Color(150, 165, 185));
            window.draw(secondLabel);

            sf::Text secondValue(formatValue(otherValues[index + 1].second), font, 11);
            secondValue.setPosition(panelX + 230.0f, rowY);
            secondValue.setFillColor(sf::Color(220, 225, 235));
            window.draw(secondValue);
        }

        window.display();

        ++fpsFrameCount;
        const sf::Time fpsElapsed = fpsClock.getElapsedTime();
        if (fpsElapsed.asSeconds() >= 1.0f) {
            currentFps = static_cast<float>(fpsFrameCount) / fpsElapsed.asSeconds();
            fpsFrameCount = 0;
            fpsClock.restart();
        }
    }
}