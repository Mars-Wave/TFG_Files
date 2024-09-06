#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ranges>
#include <set>

using namespace std;

int PRINT_TO_COUT = 0;


int findKLevelOfBranch(const std::vector<vector<int> > &vecOfVec) {
    int minNonZero = std::numeric_limits<int>::max(); // Start with the largest possible value
    for (const auto &vec: vecOfVec) {
        for (const auto &num: vec) {
            if (num != 0 && num < minNonZero) {
                minNonZero = num;
            }
        }
    }

    return minNonZero == std::numeric_limits<int>::max() ? -1 : minNonZero;
    // Return -1 if no non-zero element was found
}

int findKLevelOfNextBranch(const vector<int> &vec) {
    int minNonZero = std::numeric_limits<int>::max(); // Start with the largest possible value
    for (const auto &num: vec) {
        if (num != 0 && num < minNonZero) {
            minNonZero = num;
        }
    }

    return minNonZero == std::numeric_limits<int>::max() ? -1 : minNonZero;
    // Return -1 if no non-zero element was found
}

// Recursive function to find the next useful, valid distribution
bool findNextUsefulDistribution(const int remainingPatients, const vector<vector<int> > &currentSolution,
                                const vector<int> &roomCapacities, set<vector<vector<int> > > &exploredDistributions,
                                const int &kLevelToBeat, vector<int> &nextDistribution) {
    if (remainingPatients != 0 && remainingPatients <= kLevelToBeat) {
        return false; // Optimization: No use in exploring this one further, since split resulted in bad K-Level already
    }

    // If the current distribution sums up to the remainingPatients and is within the number of rooms and is useful (all elements > minKlevel)
    if (remainingPatients == 0 && nextDistribution.size() <= roomCapacities.size() &&
        findKLevelOfNextBranch(nextDistribution) > kLevelToBeat) {
        // Check if this distribution has not been used already
        vector<vector<int> > currentDistWrapper = currentSolution; // Match the set type for some branch node
        currentDistWrapper.push_back(nextDistribution);
        if (!exploredDistributions.contains(currentDistWrapper)) {
            return true; // Return true to indicate a valid distribution has been found
        }
        return false;
    }

    // If the number of rooms used exceeds the available rooms, return
    if (nextDistribution.size() >= roomCapacities.size()) {
        return false;
    }

    // Explore further patient distributions
    for (int i = std::min(remainingPatients, roomCapacities[nextDistribution.size()]); i >= 0; --i) {
        //only check promising patient splits
        if (i == 0 || i > kLevelToBeat) {
            // Include i patients in the current distribution
            nextDistribution.push_back(i);
            // Recursively find distributions with the remaining patients
            if (findNextUsefulDistribution(remainingPatients - i, currentSolution, roomCapacities,
                                           exploredDistributions, kLevelToBeat, nextDistribution)) {
                return true; // Return immediately if a valid distribution is found
            }
            // Backtrack: remove the last element and continue
            nextDistribution.pop_back();
        }
    }

    return false; // No valid distribution was found in this path
}

std::pair<vector<int>, bool> getNextUsefulPatientDistribution(const vector<vector<int> > &currentSolution,
                                                              set<vector<vector<int> > > &exploredDistributions,
                                                              const int kLevelToBeat, const int totalPatients,
                                                              const vector<int> &roomCapacities) {
    // Current setting is not promising and different from root
    if (findKLevelOfBranch(currentSolution) != -1 &&
        findKLevelOfBranch(currentSolution) <= kLevelToBeat) {
        return {{}, false};
    }

    vector<int> nextDistribution = {};

    // Start the recursive distribution generation, including zero patients
    bool foundNewDistribution = findNextUsefulDistribution(totalPatients, currentSolution,
                                                           roomCapacities, exploredDistributions, kLevelToBeat,
                                                           nextDistribution);

    return {nextDistribution, foundNewDistribution};
}

// Function to calculate the remaining capacities
vector<int> calculateNextCapacities(const vector<vector<int> > &currentSolution, const vector<int> &roomCapacities) {
    vector<int> remainingCapacities = roomCapacities;

    // Calculate the sum of each column in the currentSolution matrix
    for (int col = 0; col < roomCapacities.size(); ++col) {
        int columnSum = 0;

        for (const auto &row: currentSolution) {
            if (col < row.size()) { columnSum += row[col]; }
        }
        // Calculate remaining capacity for each room, invariantly >= 0 because of currentSolution's construction
        remainingCapacities[col] = roomCapacities[col] - columnSum;
    }

    return remainingCapacities;
}

void prettyPrintBestSolution(const tuple<vector<vector<int> >, int, int> &bestSolution, size_t maxColumns) {
    if (std::get<0>(bestSolution).empty()) {
        cout << "No valid patient distributions exist with the given total and room capacities." << endl;
    } else {
        cout << "Valid patient distributions for the given total and room capacities:" << endl;
        cout << "\t\t K-Level of the solution: " << std::get<1>(bestSolution) << endl;

        // Determine the number of rows (groups)
        size_t numGroups = std::get<0>(bestSolution).size();

        // Print header row
        cout << "Patient Groups / Rooms:" << endl;
        for (size_t r = 0; r < maxColumns; ++r) {
            cout << " Room " << r + 1 << " ";
        }
        cout << endl;

        // Print each patient's group distribution
        for (size_t g = 0; g < numGroups; ++g) {
            cout << "Group " << g + 1 << ":";
            for (size_t r = 0; r < maxColumns; ++r) {
                if (r < std::get<0>(bestSolution)[g].size()) {
                    // Print actual value
                    cout << " " << setw(3) << std::get<0>(bestSolution)[g][r] << " ";
                } else {
                    // Print zero for missing elements
                    cout << " " << setw(3) << 0 << " ";
                }
            }
            cout << endl;
        }
        cout << "Number of nodes visited: " << std::get<2>(bestSolution) << endl;
    }
}

void recursiveKanonymousSearch(vector<vector<int> > &bestSolution, vector<vector<int> > &currentSolution,
                               int &globalMin, const vector<int> &patientGroups,
                               const vector<int> &roomCapacities, const int maximumKLevelPossible,
                               set<vector<vector<int> > > &exploredDistributions,
                               int groupIndex = 0) {
    if (globalMin >= maximumKLevelPossible) {
        // Termination condition: K-level of assignments equals the smallest group in max bound
        // Backtrack condition: We have already found a solution of a K-level and we need to backtrack until the first promising node
        return;
    }

    if (groupIndex >= patientGroups.size()) {
        int currentMin = findKLevelOfBranch(currentSolution);
        if (currentMin > globalMin) {
            globalMin = currentMin;
            bestSolution = currentSolution;
        }
        return;
    }

    // Continue exploring all possible valid and unexplored distributions for the current group
    while (globalMin < maximumKLevelPossible) {
        // Get the next valid and unexplored distribution for the current group
        auto nextBranch = getNextUsefulPatientDistribution(
            currentSolution, exploredDistributions, globalMin, patientGroups[groupIndex],
            calculateNextCapacities(currentSolution, roomCapacities));

        // If there is no more promising and unexplored distributions, exit the loop
        if (!nextBranch.second) {
            // Mark the current distribution as explored
            exploredDistributions.emplace(currentSolution);
            break;
        }

        // Push the next valid distribution to the current solution
        currentSolution.push_back(nextBranch.first);

        // Mark the current distribution as explored
        exploredDistributions.emplace(currentSolution);

        // Recursively search for the next level
        recursiveKanonymousSearch(bestSolution, currentSolution, globalMin, patientGroups,
                                  roomCapacities, maximumKLevelPossible, exploredDistributions, groupIndex + 1);

        // Backtrack: remove the last element and continue
        currentSolution.pop_back();
    }
}

int findSmallestGroupSplitThatFits(std::vector<int> patientGroups, const std::vector<int> &roomCapacities) {
    while (true) {
        // Sort patient groups in ascending order
        ranges::sort(patientGroups);

        for (int capacity: roomCapacities) {
            if (patientGroups.front() <= capacity) {
                return patientGroups.front();
            }
        }
        std::vector<int> splitGroups;

        int half1 = patientGroups.front() / 2;
        int half2 = patientGroups.front() - half1;

        // Add the new subgroups to the splitGroups vector
        splitGroups.push_back(half1);
        splitGroups.push_back(half2);


        // Replace the original patientGroups with splitGroups
        patientGroups = splitGroups;
    }
}


tuple<vector<vector<int> >, int, int> findKanonymousDistribution(const vector<int> &patientGroups,
                                                                 const vector<int> &roomCapacities) {
    vector<vector<int> > bestSolution; //modified in recursion, return value
    int globalMin = -2; //modified in recursion
    set<vector<vector<int> > > exploredDistributions = {}; //modified in recursion
    {
        vector<vector<int> > currentSolution = {}; //modified in recursion
        recursiveKanonymousSearch(bestSolution, currentSolution, globalMin, patientGroups, roomCapacities,
                                  findSmallestGroupSplitThatFits(patientGroups, roomCapacities), exploredDistributions);
    }

    if (PRINT_TO_COUT == 1) {
        if (findSmallestGroupSplitThatFits(patientGroups, roomCapacities) == globalMin) {
            cout << "------- TERMINATION CONDITION MET ----------" << endl << endl;
        } else {
            cout << "\n----- Termination condition was too relaxed, at : " <<
                    findSmallestGroupSplitThatFits(patientGroups, roomCapacities) << "----------" << endl << endl;
        }
    }

    return {bestSolution, globalMin, exploredDistributions.size()};
}


int parseInput(char *argv[], vector<int> &groups, vector<int> &rooms) {
    // Arguments are provided in the following order:
    // MAXIMUM_POSSIBLE_KLEVEL
    // NUMBER_ROOMS
    // PATIENTS_GROUPS
    // room capacities...
    // patient groups...

    int maxPossibleKLevel = 0;

    // if argv[2] == 0 we assume we are in batch mode and printing will slow the process down
    PRINT_TO_COUT = argv[2][0] - '0';
    int numRooms, numGroups;
    // Open the input file
    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error: Could not open input file.");
    }

    // Parse values from the file
    inputFile >> maxPossibleKLevel >> numRooms >>
            numGroups;

    // Parse room capacities
    rooms.resize(numRooms);
    for (int i = 0; i < numRooms; ++i) {
        inputFile >> rooms[i];
    }

    // Parse patients
    groups.resize(numGroups);
    for (int i = 0; i < numGroups; ++i) {
        inputFile >> groups[i];
    }

    inputFile.close();
    return maxPossibleKLevel;
}


void appendDataToCSV(const std::string &filename, const std::vector<std::string> &newData) {
    std::ofstream file(filename,
                       std::ios_base::app); // Open file in append mode

    if (!file.is_open()) {
        std::cerr << "Error: Could not open " + filename + " for appending."
                << std::endl;
        return;
    }
    for (const auto &data: newData) {
        file << data << ",";
    }
    file << std::endl;
    file.close();
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: <filename> <printing>" << std::endl;
        return 1;
    }
    vector<int> patientGroups;
    vector<int> roomCapacities;

    int maxPossibleKlevel = parseInput(argv, patientGroups, roomCapacities);

    int totalPatients = accumulate(patientGroups.begin(), patientGroups.end(), 0);
    int totalCapacity = accumulate(roomCapacities.begin(), roomCapacities.end(), 0);

    // Validate if total patients can be accommodated
    if (totalPatients > totalCapacity) {
        cout << "Error: Total number of patients exceeds total room capacities." << endl;
        return -1;
    }

    // Start timer
    auto startTime = std::chrono::steady_clock::now();

    auto bestSolution = findKanonymousDistribution(patientGroups, roomCapacities);
    // Stop timer
    auto endTime = std::chrono::steady_clock::now();

    // Calculate the elapsed time
    auto elapsedTime =
            std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime)
            .count();

    if (PRINT_TO_COUT) {
        prettyPrintBestSolution(bestSolution, roomCapacities.size());
    } else {
        std::string input = argv[1];
        std::string baseName = input.substr(0, input.find_last_of('.'));
        std::string output = "results" + baseName + "EXHAUST.csv";
        appendDataToCSV(output, {
                            std::to_string(maxPossibleKlevel),
                            std::to_string(std::get<1>(bestSolution)),
                            std::to_string(std::get<2>(bestSolution)),
                            std::to_string(elapsedTime)
                        });
    }
    return 0;
}
