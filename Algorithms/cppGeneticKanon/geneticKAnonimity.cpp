#include <algorithm>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <crtdbg.h>
#include <unordered_set>

// Global variables

// Best K-Level possible is given as an argument
int MAXIMUM_POSSIBLE_KLEVEL;

// Define the static capacity of a room (e.g., maximum capacity).
int ROOM_CAPACITY;

// Define the group range of patients.
int PATIENTS_GROUPS;

// Define the total number of rooms available (e.g., for allocation).
int NUMBER_ROOMS;

// Define the maximum number of generations
int MAX_GENERATIONS;

// Define the population size
int POPULATION_SIZE;

// Define the chance of mutating of a gene
int MUTATION_CHANCE;

// Allow printing lines in stdout 0=no, 1=yes
int PRINT_TO_COUT;

// Result mode
// 0=mean of 10 points for 1 result of total runtime + K-level of best + Overall violations
// 1=growth evolution of best in every generation with K-Level and violations, with generations until convergence
int RESULT_MODE;

// MetaHeuristicVersion or Random version 0=Random, 1=Meta
int SOLUTION_TYPE;

//output file
std::string OUTPUT;

// Aliases for code clarity
enum class SolutionType : int {
    RandomType = 0,
    HeuristicGuidedType = 1,
    SelfRepairingHeuristicGuidedType = 2
};

// Data definitions
struct PatientGroup {
    int amountOfPatients;
};

struct Room {
    int capacity;
};

//Problem structures

struct FitnessScore {
    int roomsOverMaxCapacity, groupsOverMaxCapacity, levelOfK;

    //Deterministic Fitness Scoring

    // Define the less than operator to compare FitnessScore objects
    // i.e: a < b implies a is less fit than b
    bool operator<(const FitnessScore &other) const {
        // Compare based on roomsOverMaxCapacity
        if (roomsOverMaxCapacity != other.roomsOverMaxCapacity) {
            return roomsOverMaxCapacity > other.roomsOverMaxCapacity;
        }
        // Compare based on groupsOverMaxCapacity
        if (groupsOverMaxCapacity != other.groupsOverMaxCapacity) {
            return groupsOverMaxCapacity > other.groupsOverMaxCapacity;
        }
        // If groupsOverMaxCapacity is equal, compare based on levelOfK
        return levelOfK < other.levelOfK;
    }
};

struct ProblemContext {
    static std::vector<Room> roomCapacities;
    static std::vector<PatientGroup> patientGroups;
};

// Initialize static members
std::vector<Room> ProblemContext::roomCapacities;
std::vector<PatientGroup> ProblemContext::patientGroups;

class Solution {
protected:
    std::vector<std::vector<int> > assignments;

    std::vector<int> peoplePerRoom;

    std::vector<int> patientsPerGroup;

    std::priority_queue<std::tuple<int, int, int>,
        std::vector<std::tuple<int, int, int> >,
        std::greater<> > kLevelminHeap;

    int minPeople;
    int roomsOverCapacityCount;
    int wrongGroupsCount;


    /// Initialize the min-heap with non-zero assignment values
    void initKLevelMinHeap() {
        for (size_t i = 0; i < assignments.size(); ++i) {
            for (size_t j = 0; j < assignments[i].size(); ++j) {
                if (assignments[i][j] > 0) {
                    kLevelminHeap.emplace(assignments[i][j], i, j);
                }
            }
        }
    }

    void clearMinHeap() {
        std::priority_queue<std::tuple<int, int, int>,
            std::vector<std::tuple<int, int, int> >,
            std::greater<> > newHeap;

        while (!kLevelminHeap.empty()) {
            auto [value, groupId, roomId] = kLevelminHeap.top();
            if (assignments[groupId][roomId] == value) {
                newHeap.push(kLevelminHeap.top());
            }
            kLevelminHeap.pop();
        }

        kLevelminHeap = std::move(newHeap);
    }

    void updatePeoplePerRoom(int roomId, int delta) {
        int oldCount = peoplePerRoom[roomId];
        peoplePerRoom[roomId] += delta;
        if (oldCount <= getRoomCapacity(roomId) &&
            peoplePerRoom[roomId] > getRoomCapacity(roomId)) {
            // We exceeded capacity in this room
            roomsOverCapacityCount++;
        } else if (oldCount > getRoomCapacity(roomId) &&
                   peoplePerRoom[roomId] <= getRoomCapacity(roomId)) {
            // We went from exceeding capacity to being in range
            roomsOverCapacityCount--;
        }
    }

    void updatePatientsPerGroup(int groupId, int delta) {
        int oldGroupCount = patientsPerGroup[groupId];
        patientsPerGroup[groupId] += delta;
        if (oldGroupCount != getPatientGroupCount(groupId) &&
            patientsPerGroup[groupId] == getPatientGroupCount(groupId)) {
            wrongGroupsCount--;
        } else if (oldGroupCount == getPatientGroupCount(groupId) &&
                   patientsPerGroup[groupId] != getPatientGroupCount(groupId)) {
            wrongGroupsCount++;
        }
    }

    static int getRoomCapacity(int roomId) {
        return ProblemContext::roomCapacities[roomId].capacity;
    }

    static int getPatientGroupCount(int groupId) {
        return ProblemContext::patientGroups[groupId].amountOfPatients;
    }

    void updateKLevel() {
        // Clear min entries from the heap
        while (!kLevelminHeap.empty()) {
            auto [value, row, col] = kLevelminHeap.top();
            if (assignments[row][col] == value) {
                minPeople = value;
                return;
            }
            kLevelminHeap.pop();
        }
    }

    void initAmountInGroupToRoom(int groupId, int roomId, int amount = -1) {
        if (amount == -1) {
            amount =
                    rand() %
                    (std::min(getRoomCapacity(roomId),
                              getPatientGroupCount(groupId)) + 1);
        }
        if (amount > 0) {
            // Update K-level
            kLevelminHeap.emplace(amount, groupId, roomId);
            updatePeoplePerRoom(roomId, amount);
            updatePatientsPerGroup(groupId, amount);
        }
        // Update the assignments matrix
        assignments[groupId][roomId] = amount;
    }


    void updateDataStructures(int groupId, int roomId, int amount) {
        int oldAmount = assignments[groupId][roomId];
        assignments[groupId][roomId] = amount;
        if (amount > 0) {
            kLevelminHeap.emplace(amount, groupId, roomId);
            // Update K-level
            updateKLevel();
        }
        updatePeoplePerRoom(roomId, amount - oldAmount);
        updatePatientsPerGroup(groupId, amount - oldAmount);
    }

public:
    // Initialize the solution with rooms and groups
    Solution()
        : assignments(PATIENTS_GROUPS, std::vector<int>(NUMBER_ROOMS, 0)),
          peoplePerRoom(NUMBER_ROOMS, 0), patientsPerGroup(PATIENTS_GROUPS, 0),
          minPeople(INT_MAX), roomsOverCapacityCount(0), wrongGroupsCount(PATIENTS_GROUPS), kLevelminHeap(
              std::priority_queue<std::tuple<int, int, int>,
                  std::vector<std::tuple<int, int, int> >,
                  std::greater<> >{}) {
    }

    virtual ~Solution() = default;

    void initializeSolution() {
        // Initialize internal data
        for (int groupId = 0; groupId < PATIENTS_GROUPS; ++groupId) {
            for (int roomId = 0; roomId < NUMBER_ROOMS; ++roomId) {
                initAmountInGroupToRoom(groupId, roomId); // Random init
            }
        }

        // Update KLevel after initialization
        updateKLevel();
    }

    void assignAmountInGroupToRoom(int groupId, int roomId, int amount) {
        updateDataStructures(groupId, roomId, amount);
    }

    virtual void repairSomeGroupInGene() {
    };

    virtual void mutateAmountInGroupAndRoom(int groupId, int roomId) = 0;

    ///auxiliary function for crossover
  [[nodiscard]]  virtual Solution *clone() const = 0;

[[nodiscard]]    int getPeopleFromGroupInRoom(int group, int room) const {
        return assignments[group][room];
    }

    /// Get the number of people in a specific room
  [[nodiscard]]  int getPeopleInRoom(int roomId) const { return peoplePerRoom[roomId]; }

    /// Get the number of rooms that are over capacity.
  [[nodiscard]]  int getRoomsOverCapacityCount() const { return roomsOverCapacityCount; }

    /// Get the number of groups that have an incorrect number of patients.
   [[nodiscard]] int getWrongGroupsCount() const { return wrongGroupsCount; }

    /// Get the minimum number of people from a group assigned to any
    /// room (K-level).
    [[nodiscard]] int getKLevel() const { return minPeople; }

    /// Get defined fitness score
    [[nodiscard]] FitnessScore getFitness() const {
        return {roomsOverCapacityCount, wrongGroupsCount, getKLevel()};
    }

    //printer
    void printSolution() const {
        // Print assignments matrix
        std::cout << "Assignments Matrix:" << std::endl;
        std::cout << "-------------------" << std::endl;
        for (int groupId = 0; groupId < PATIENTS_GROUPS; ++groupId) {
            std::cout << "Group " << std::setw(3) << groupId + 1 << ": ";
            for (int roomId = 0; roomId < NUMBER_ROOMS; ++roomId) {
                std::cout << std::setw(3) << assignments[groupId][roomId] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        // Print other statistics
        std::cout << "Statistics:" << std::endl;
        std::cout << "-----------" << std::endl;
        std::cout << "Bad Rooms: " << roomsOverCapacityCount << std::endl;
        std::cout << "Bad Groups: " << wrongGroupsCount << std::endl;
        std::cout << "K-Level: " << minPeople << std::endl;
        std::cout << std::endl;
    }
};

/** Override mutateAmountInGroupAndRoom in a standard, random approach
     * In this approach there is:
     *      - a 100% probability to make the assignation a random split of the group
     */
class RandomSolution : public Solution {
public:
    using Solution::Solution; // Inherit constructors if needed

    [[nodiscard]] Solution *clone() const override {
        return new RandomSolution(*this);
    }

    //no repair
    void repairSomeGroupInGene(int groupId, int roomId) {
        clearMinHeap();
    }

    // Override mutateAmountInGroupAndRoom for random approach
    void mutateAmountInGroupAndRoom(int groupId, int roomId) override {
        int amount =
                std::rand() %
                (std::min(getRoomCapacity(roomId),
                          getPatientGroupCount(groupId)) + 1);

        assignAmountInGroupToRoom(groupId, roomId, amount);
    }
};

/** Overrides mutateAmountInGroupAndRoom for HeuristicGuidedType approach
     * In this approach there is:
     *      - a 25% probability to make the assignation 0
     *      - a 50% probability to make the assignation the max group count for that group
     *      - a 25% probability to make the assignation a random split of the group
     */
class HeuristicGuidedSolution : public Solution {
public:
    using Solution::Solution; // Inherit constructors if needed

    [[nodiscard]] Solution *clone() const override {
        return new HeuristicGuidedSolution(*this);
    }

    //no repair
    void repairSomeGroupInGene(int groupId, int roomId) {
        clearMinHeap();
    }

    void mutateAmountInGroupAndRoom(int groupId, int roomId) override {
        // Generate a random number between 0 and 99
        int randValue = std::rand() % 100;

        int amount;
        if (randValue < 25) {
            // 25% probability
            amount = 0;
        } else if (randValue < 75) {
            // 50% probability (25% + 50%)
            amount = getPatientGroupCount(groupId);
        } else {
            // 25% probability
            amount = std::rand() % (std::min(getRoomCapacity(roomId), getPatientGroupCount(groupId)) + 1);
        }

        assignAmountInGroupToRoom(groupId, roomId, amount);
    }
};


/** Overrides mutateAmountInGroupAndRoom for SelfRepairingHeuristicGuidedSolution approach
     * In this approach there is:
     *      - a 25% probability to make the assignation 0
     *      - a 50% probability to make the assignation the max group count for that group
     *      - a 25% probability to make the assignation a random split of the group
     *
     * And a reparation of the matrix ensuring that, after a mutation, there never is:
     *      - a group assigned unequal to the statement
     *      - a room exceeding capacity
     */
class SelfRepairingHeuristicGuidedSolution : public Solution {
private:
    int getTotalGroupDifferenceForAssignment(int groupId, int roomId) {
        int groupCount = getPatientGroupCount(groupId);
        int currentGroupSum = patientsPerGroup[groupId];
        if (currentGroupSum == groupCount) {
            return -1;
        }

        return std::max(std::min(std::max(0, assignments[groupId][roomId] - (currentGroupSum - groupCount)),
                                 getRoomCapacity(roomId) - peoplePerRoom[roomId]), 0);
    }

    int getRoomDecrementForGroup(int groupId, int roomId) {
        int roomCapacity = getRoomCapacity(roomId);
        int currentRoomSum = peoplePerRoom[roomId];
        if (currentRoomSum <= roomCapacity) {
            return -1;
        }

        return std::rand() % (std::max(0, assignments[groupId][roomId] - (currentRoomSum - roomCapacity)) +
                              1);
    }

public:
    using Solution::Solution; // Inherit constructors if needed

    void repairSomeGroupInGene() override {
        std::vector<int> badGroups;
        for (int i = 0; i < PATIENTS_GROUPS; i++) {
            if (getPatientGroupCount(i) != patientsPerGroup
                [i]) {
                badGroups.push_back(i);
            }
        }
        if (badGroups.empty()) {
            clearMinHeap();
            return;
        }
        int wrongGroup = badGroups[rand() % badGroups.size()];

        if (getPatientGroupCount(wrongGroup) > patientsPerGroup[wrongGroup]) {
            // Underassignment case
            // Collect rooms with available capacity
            std::vector<int> availableRooms;
            for (int r = 0; r < NUMBER_ROOMS; r++) {
                if (getRoomCapacity(r) > peoplePerRoom[r]) {
                    availableRooms.push_back(r);
                }
            }
            while (getPatientGroupCount(wrongGroup) != patientsPerGroup[wrongGroup]) {
                if (availableRooms.empty()) {
                    // No available rooms to assign patients to, exit the loop
                    break;
                }

                int targetRoom = availableRooms[rand() % availableRooms.size()];
                int availableCapacity = std::max(getRoomCapacity(targetRoom) - peoplePerRoom[targetRoom], 0);

                // Determine how many patients we can move to the target room
                int patientsToMove = std::min(getPatientGroupCount(wrongGroup) - patientsPerGroup[wrongGroup],
                                              availableCapacity);

                // Update assignments and data structures
                updateDataStructures(wrongGroup, targetRoom, assignments[wrongGroup][targetRoom] + patientsToMove);

                // If the room is now full, remove it from available rooms
                if (peoplePerRoom[targetRoom] >= getRoomCapacity(targetRoom)) {
                    availableRooms.erase(std::remove(availableRooms.begin(), availableRooms.end(), targetRoom),
                                         availableRooms.end());
                }
            }
        } else if (getPatientGroupCount(wrongGroup) < patientsPerGroup[wrongGroup]) {
            // Overassignment case
            // Collect rooms with patients on it
            std::vector<int> roomsOccupied;
            for (int r = 0; r < NUMBER_ROOMS; r++) {
                for (int g = 0; g < PATIENTS_GROUPS; g++) {
                    if (assignments[g][r] > 0) {
                        roomsOccupied.push_back(r);
                        break;
                    }
                }
            }
            while (getPatientGroupCount(wrongGroup) != patientsPerGroup[wrongGroup]) {
                if (roomsOccupied.empty()) {
                    // No available rooms to subtract patients from, exit the loop
                    break;
                }

                int targetRoom = roomsOccupied[rand() % roomsOccupied.size()];

                // Determine how many patients we can erase in the target room
                int newAmount = std::max(
                    assignments[wrongGroup][targetRoom] - (
                        patientsPerGroup[wrongGroup] - getPatientGroupCount(wrongGroup)), 0);

                // Update assignments and data structures
                updateDataStructures(wrongGroup, targetRoom, newAmount);

                // If the room is now full, remove it from available rooms
                if (assignments[wrongGroup][targetRoom] == 0) {
                    roomsOccupied.erase(std::remove(roomsOccupied.begin(), roomsOccupied.end(), targetRoom),
                                        roomsOccupied.end());
                }
            }
        }
        clearMinHeap();
    }

    [[nodiscard]] Solution *clone() const override {
        return new SelfRepairingHeuristicGuidedSolution(*this);
    }

    void mutateAmountInGroupAndRoom(int groupId, int roomId) override {
        // Generate a random number between 0 and 99
        int randValue = std::rand() % 100;

        int amount;
        if (randValue < 25) {
            // 25% probability
            amount = 0;
        } else if (randValue < 75) {
            // 50% probability (25% + 50%)
            amount = getPatientGroupCount(groupId);
        } else {
            // 25% probability
            amount = std::rand() % (std::min(getRoomCapacity(roomId), getPatientGroupCount(groupId)) + 1);
        }

        assignAmountInGroupToRoom(groupId, roomId, amount);
    }
};

[[nodiscard]] Solution *SelectParent(std::vector<Solution *> &population) {
    // Implementation of tournament selection

    // Defintion of the size of the tournament (2 for a binary tournament).
    const int tournamentSize = 2;

    // Randomly select individuals for the tournament.
    auto bestParent = population[rand() %
                                 population.size()];
    FitnessScore bestFitness = bestParent->getFitness();

    for (int i = 0; i < tournamentSize - 1; ++i) {
        int randomIndex =
                rand() %
                population.size(); // Randomly choose an individual from the population.
        const FitnessScore currentFitness = population[randomIndex]->getFitness();

        if (bestFitness < currentFitness) {
            bestParent = population[randomIndex];
            bestFitness = currentFitness;
        }
    }

    return bestParent;
}

[[nodiscard]] Solution *Crossover(const Solution * &parent1, const Solution * &parent2) {
    // Implementation of crossover method, two-point
    // crossover, to create a child solution from the parents. Returns the child
    // solution as a vector of Assignments
    Solution *child = parent1->clone();


    // Randomly choose two distinct columns for the crossover.
    int point1 = rand() % NUMBER_ROOMS;
    int point2 = rand() % NUMBER_ROOMS;

    // Ensure point1 is less than point2.
    if (point1 > point2) {
        std::swap(point1, point2);
    }

    // Replace groups from parent1 with groups from point1 to point2 from parent2
    for (int j = point1; j < std::min(point2, NUMBER_ROOMS); j++) {
        for (int i = 0; i < PATIENTS_GROUPS; i++) {
            int amount = parent2->getPeopleFromGroupInRoom(i, j);
            child->assignAmountInGroupToRoom(i, j, amount);
        }
    }

    return child;
}

void Mutate(Solution *child) {
    // Implement a mutation method to introduce random changes to the child
    // solution.
    // Randomly choose an assignment to mutate.
    if (rand() % 100 <= MUTATION_CHANCE) {
        int mutationI = rand() % PATIENTS_GROUPS;
        int mutationJ = rand() % NUMBER_ROOMS;

        child->mutateAmountInGroupAndRoom(mutationI, mutationJ);
    }
}

[[nodiscard]] Solution *GetBestSolution(std::vector<Solution *> &population) {
    // Find and return the best solution from the population based on fitness.
    // The best solution is the one with the best fitness score.

    Solution *bestSolution = population[0];
    FitnessScore bestFitness = bestSolution->getFitness();

    for (Solution * &solution: population) {
        FitnessScore fitness = solution->getFitness();

        if (bestFitness < fitness) {
            bestFitness = fitness;
            bestSolution = solution;
        }
    }
    return bestSolution;
}

void ReplaceWeakest(std::vector<Solution *> &population, Solution * &child) {
    // Compare the fitness of the child solution with the fitness of the weakest
    // solution in the population. If the child's fitness is better, replace the
    // weakest solution with the child solution.

    // Evaluate fitness
    FitnessScore childFitness = child->getFitness();

    // Find the index of the weakest solution in the population.
    FitnessScore weakestFitness = population[0]->getFitness();
    int weakestIndex = 0;

    for (int i = 1; i < population.size(); ++i) {
        FitnessScore fitness = population[i]->getFitness();
        if (fitness < weakestFitness) {
            weakestFitness = fitness;
            weakestIndex = i;
        }
    }

    // Compare the child's fitness with the weakest solution's fitness.
    if (weakestFitness < childFitness) {
        // Replace the weakest solution with the child solution.
        delete population[weakestIndex];
        population[weakestIndex] = child;
    }
}

/// Dump results in a new line of results.csv
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

[[nodiscard]] Solution *GeneticAlgorithm(std::vector<Solution *> &population) {
    // Main Genetic Algorithm Loop
    for (int generation = 0; generation < MAX_GENERATIONS; ++generation) {
        // Selection
        const auto *parent1 = SelectParent(population);
        const auto *parent2 = SelectParent(population);

        // Crossover
        Solution *child = Crossover(parent1, parent2);

        // Mutation
        Mutate(child);

        //Gene repair
        child->repairSomeGroupInGene();

        // Update population
        ReplaceWeakest(population, child);

        if (RESULT_MODE == 1 && (generation % 10 == 0)) {
            const auto best = GetBestSolution(population);
            appendDataToCSV(OUTPUT, {
                                std::to_string(generation),
                                std::to_string(MAXIMUM_POSSIBLE_KLEVEL),
                                std::to_string(best->getKLevel()),
                                std::to_string(best->getRoomsOverCapacityCount()),
                                std::to_string(best->getWrongGroupsCount())
                            });
        }
    }
    // Output the best solution
    return GetBestSolution(population);
}

void PrintGeneticAlgorithmContext() {
    std::cout << "Context of the Genetic Algorithm:" << std::endl;
    std::cout << "\tMAX_GENERATIONS: " << MAX_GENERATIONS
            << " (Maximum number of generations)" << std::endl;
    std::cout << "\tPOPULATION_SIZE: " << POPULATION_SIZE << " (Population size)"
            << std::endl;
    std::cout << "\tMUTATION_CHANCE: " << MUTATION_CHANCE << "% (Chance of mutation 1-100)" << std::endl;
}

void PrintGeneticAlorithmVariables() {
    std::cout << "Variables of the Genetic Algorithm:" << std::endl;

    std::cout << "\tNUMBER_ROOMS: " << NUMBER_ROOMS
            << " (Total number of available rooms for allocation)" << std::endl;
    for (int i = 0; i < NUMBER_ROOMS; ++i) {
        std::cout << "\t\tCapacity for room " << i + 1 << ": " << ProblemContext::roomCapacities[i].capacity <<
                std::endl;
    }

    std::cout << "\tPATIENTS_GROUPS: " << PATIENTS_GROUPS
            << " (Group range of patients)" << std::endl;
    for (int i = 0; i < PATIENTS_GROUPS; ++i) {
        std::cout << "\t\tAmount of patients for k-group " << i + 1 << ": " << ProblemContext::patientGroups[i].
                amountOfPatients << std::endl;
    }
}

void parseInput(char *argv[]) {
    // Arguments are provided in the following order:
    //  MAX_GENERATIONS
    //  POPULATION_SIZE
    // MUTATION_CHANCE
    // MAXIMUM_POSSIBLE_KLEVEL
    // NUMBER_ROOMS
    // PATIENTS_GROUPS
    // room capacities...
    // patient groups...


    // if argv[2] == 0 we assume we are in batch mode and printing will slow the process down
    PRINT_TO_COUT = argv[2][0] - '0';

    // argv[3] is solution version
    SOLUTION_TYPE = argv[3][0] - '0';

    // argv[4] == is result mode to csv
    RESULT_MODE = argv[4][0] - '0';

    // Open the input file
    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error: Could not open input file.");
    }

    // Parse values from the file
    inputFile >> MAX_GENERATIONS >> POPULATION_SIZE >> MUTATION_CHANCE >> MAXIMUM_POSSIBLE_KLEVEL >> NUMBER_ROOMS >>
            PATIENTS_GROUPS;

    if (PRINT_TO_COUT) {
        PrintGeneticAlgorithmContext();
    }

    // Parse room capacities
    ProblemContext::roomCapacities.resize(NUMBER_ROOMS);
    for (int i = 0; i < NUMBER_ROOMS; ++i) {
        inputFile >> ProblemContext::roomCapacities[i].capacity;
    }

    // Parse patients
    ProblemContext::patientGroups.resize(PATIENTS_GROUPS);
    for (int i = 0; i < PATIENTS_GROUPS; ++i) {
        inputFile >> ProblemContext::patientGroups[i].amountOfPatients;
    }

    if (PRINT_TO_COUT) {
        PrintGeneticAlorithmVariables();
    }

    inputFile.close();
}

int main(int argc, char *argv[]) {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    // Parse input file if provided for results.csv populating
    if (argc > 4) {
        try {
            parseInput(argv);
        } catch (const std::runtime_error &e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    } else {
        std::cout << "Usage: <filename> <printing> <solutiontype> <resultmode>" << std::endl;
        return 1;
    }

    std::string input = argv[1];
    std::string baseName = input.substr(0, input.find_last_of('.'));
    OUTPUT = "results" + baseName + "_V" + std::to_string(SOLUTION_TYPE) + (RESULT_MODE == 1 ? "Evo" : "mean10")
             + ".csv";


    if (RESULT_MODE == 0) {
        long double totalTime = 0;
        double cumulKLevel = 0;
        double cumulViolations = 0;

        std::vector<Solution *> population(POPULATION_SIZE, nullptr);

        for (int i = 0; i < 10; i++) {
            for (size_t i = 0; i < POPULATION_SIZE; ++i) {
                if (SOLUTION_TYPE == static_cast<int>(SolutionType::RandomType)) {
                    population[i] = new RandomSolution();
                } else if (SOLUTION_TYPE == static_cast<int>(SolutionType::HeuristicGuidedType)) {
                    population[i] = new HeuristicGuidedSolution();
                } else if (SOLUTION_TYPE == static_cast<int>(
                               SolutionType::SelfRepairingHeuristicGuidedType)) {
                    population[i] = new SelfRepairingHeuristicGuidedSolution();
                } else {
                    std::cerr << "No solution type for argument provided" << std::endl;
                    return 1;
                }
            }

            for (const auto &sol: population) {
                sol->initializeSolution();
            }

            // Start timer
            auto startTime = std::chrono::steady_clock::now();

            // Run the genetic algorithm
            const Solution *bestSolution = GeneticAlgorithm(population);

            // Stop timer
            auto endTime = std::chrono::steady_clock::now();

            // Calculate the elapsed time
            const auto elapsedTime =
                    std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime)
                    .count();

            totalTime += elapsedTime;
            cumulKLevel += bestSolution->getKLevel();
            cumulViolations += bestSolution->getWrongGroupsCount() + bestSolution->getRoomsOverCapacityCount();

            if (PRINT_TO_COUT) {
                std::cout << "Best solution iter #" << std::to_string(i) << ":" << std::endl << "-------------------" <<
                        std::endl;
                bestSolution->printSolution();
            }

            // Clean up allocated memory
            for (auto &sol: population) {
                delete sol;
                sol = nullptr; // Avoid dangling pointers
            }
        }
        appendDataToCSV(OUTPUT, {
                            std::to_string(MAXIMUM_POSSIBLE_KLEVEL),
                            std::to_string(cumulKLevel / 10),
                            std::to_string(cumulViolations / 10),
                            std::to_string(totalTime / 10)
                        });
    }
    if (RESULT_MODE == 1) {
        std::vector<Solution *> population;

        for (size_t i = 0; i < POPULATION_SIZE; ++i) {
            if (SOLUTION_TYPE == static_cast<int>(SolutionType::RandomType)) {
                population.push_back(new RandomSolution());
            } else if (SOLUTION_TYPE == static_cast<int>(SolutionType::HeuristicGuidedType)) {
                population.push_back(new HeuristicGuidedSolution());
            } else if (SOLUTION_TYPE == static_cast<int>(
                           SolutionType::SelfRepairingHeuristicGuidedType)) {
                population.push_back(new SelfRepairingHeuristicGuidedSolution());
            } else {
                std::cerr << "No solution type for argument provided" << std::endl;
                return 1;
            }
        }

        for (auto &sol: population) {
            sol->initializeSolution();
        }

        // Run the genetic algorithm
        Solution *bestSolution = GeneticAlgorithm(population);

        if (PRINT_TO_COUT) {
            std::cout << "Best solution:" << std::endl << "--------------" <<
                    std::endl;
            bestSolution->printSolution();
        }
        // Clean up allocated memory
        for (auto solution: population) {
            delete solution;
        }
    }

    return 0;
}
