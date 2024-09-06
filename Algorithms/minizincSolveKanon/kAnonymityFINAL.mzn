include "globals.mzn";

% Data
int: kGroups;
constraint assert(kGroups > 0, "Patients are grouped");

array[1..kGroups] of int: kGroupVector;
constraint assert(forall(i in 1..kGroups)(kGroupVector[i] > 0), "Groups are positive");

int: rooms;
constraint assert(rooms > 0, "We have screening rooms");

array[1..rooms] of int: roomCapacity;
constraint assert(forall(i in 1..rooms)(roomCapacity[i] > 0), "Room Capacity is not negative");
constraint assert(sum(kGroupVector) <= sum(roomCapacity), "We can allocate all patients in our rooms");

% Decision variables
var 1..min(kGroupVector): kLevel;
array[1..kGroups, 1..rooms] of var 0..max(kGroupVector): groupRoomCounts;

% Constraints

% Ensure groupRoomCounts respects group sizes
constraint forall(g in 1..kGroups)(
  sum(r in 1..rooms)(groupRoomCounts[g, r]) == kGroupVector[g]
);

% Ensure groupRoomCounts respects room capacities
constraint forall(r in 1..rooms)(
  sum(g in 1..kGroups)(groupRoomCounts[g, r]) <= roomCapacity[r]
);

% Calculate the K-Level of the solution as min(groupRoomCounts[g,r]) non-zero for all g in groups and r in rooms
constraint kLevel == min([groupRoomCounts[g, r] | g in 1..kGroups, r in 1..rooms where groupRoomCounts[g, r] > 0]);

% Objective function to maximize the K-Level of the solution
solve maximize kLevel;

% Output the results
output [
    "K-Level of solution: \(kLevel)\n",
] ++
[ "K-Groups: \(kGroupVector)\n"
]++
[ "Room Capacities: \(roomCapacity) \n"
]++
["Group Room Counts:\n"
]++
[
  "Room \(r):\t" ++ 
  concat([show(groupRoomCounts[g, r]) ++ "\t" | g in 1..kGroups]) ++ 
  "\n" 
  | r in 1..rooms
];