#include "../include/heston.hpp"
#include "../include/csv_writer.hpp"
#include <iostream>
#include <chrono>

int main() {
    HestonParams p{};

    HestonSimulator sim(p);

    auto start = std::chrono::steady_clock::now();

    auto [S_paths, v_paths] = sim.simulate(10000, 252 * 2);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Simulate took "
              << std::chrono::duration<double>(end - start).count() << " s\n";

    write_paths_csv("../../results/heston_paths.csv", ...)

    return 0;
}