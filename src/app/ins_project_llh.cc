#include "core/navigation/ins_navigation.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    if (argc!=4) { std::cerr << "Usage: ins_project_llh origin_lat origin_lon origin_h < llh.txt\n"; return 2; }
    try {
        lightning::ins::GeoReference ref(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]));
        std::string line;
        while (std::getline(std::cin,line)) {
            if (line.empty() || line.front()=='#') continue;
            double lat,lon,h; std::string extra; std::istringstream row(line);
            if (!(row>>lat>>lon>>h) || (row>>extra)) throw std::invalid_argument("expected latitude[deg] longitude[deg] ellipsoid_height[m]");
            const auto p=ref.Forward(lat,lon,h);
            std::cout << std::setprecision(16) << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
        }
        return std::cin.bad() || !std::cout ? 1 : 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
