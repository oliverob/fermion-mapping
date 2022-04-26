
from typing import List
import numpy as np
import math


class Qubit():
    def __init__(self, parent, root_enumeration) -> None:
        self.parent = parent
        self.children = []
        self.root_enumeration = root_enumeration
        self.tree_enumeration = -1

    def add_child(self, child):
        self.children.append(child)

    def get_children(self):
        return self.children

    def get_parent(self):
        return self.parent
    
    def add_child(self, child):
        self.children.append(child)
    
    def set_tree_enumeration(self, tree_enumeration):
        self.tree_enumeration = tree_enumeration

    def get_update_qubits(self):
        if self.parent == None:
            return []
        update = self.get_parent().get_update_qubits()
        if update == []:
            return [self.parent]
        update.append(self.parent)
        return update

    def get_parity_qubits(self, i):
        parity = [child for child in self.get_children() if child.tree_enumeration < i]
        if self.get_parent() == None:
            return parity
        parent_parity = self.get_parent().get_parity_qubits(i)
        if parent_parity == []:
            return parity
        parity.extend(parent_parity)
        return parity

class IncorrectInitialisation(Exception):
    pass


class Lattice():
    def __init__(self, height, width, qubits, enumeration_x_coordinates, enumeration_y_coordinates) -> None:
        if(width*height == len(qubits) and len(qubits) == len(enumeration_x_coordinates) and len(enumeration_x_coordinates) == len(enumeration_y_coordinates)):
            self.width = width
            self.height = height
            self.lattice_array = [[0]*width for y in range(0,height)]
            self.coordinates = [[enumeration_x_coordinates[i], enumeration_y_coordinates[i]]
                                for i in range(0, len(enumeration_x_coordinates))]
            for i in range(len(qubits)):
                self.lattice_array[enumeration_y_coordinates[i]][enumeration_x_coordinates[i]] = qubits[i]
        else:
            print("Incorrect Initialisation")
            raise IncorrectInitialisation

    def get_neighbours(self, x, y) -> List[Qubit]:
        neighbours = []
        if(x > 0):
            neighbours.append([x-1,y])
        if(x < self.width - 1):
            neighbours.append([x+1,y])
        if(y > 0):
            neighbours.append([x,y-1])
        if(y < self.height - 1):
            neighbours.append([x,y+1])
        return neighbours
    
    def get_qubit_by_coordinates(self,coordinates) -> Qubit:
        x,y = coordinates
        return self.lattice_array[y][x]

    def __str__(self) -> str:
        root_enumeration_array = []
        for row in self.lattice_array:
            root_enumeration_array.append([qubit.root_enumeration for qubit in row])
        return str(root_enumeration_array)

def gates_needed_to_create_annhilate_fermionic_mode(lattice: Lattice, coordinates, fenwick_trees):
    target_qubit = lattice.get_qubit_by_coordinates(coordinates)
    print("Target qubit: ", target_qubit.root_enumeration, ", ",target_qubit.tree_enumeration)
    update_qubits = target_qubit.get_update_qubits()
    parity_qubits = target_qubit.get_parity_qubits(target_qubit.tree_enumeration)
    parity_qubits.extend([fenwick_trees[i][0] for i in range(0,target_qubit.root_enumeration)])
    return update_qubits, target_qubit, parity_qubits

def find_nonoverlapping_qubits(creation_qubit_string, annhiliation_qubit_string):
    nonoverlapping_qubits = []
    for qubit in creation_qubit_string:
        if qubit not in annhiliation_qubit_string:
            nonoverlapping_qubits.append(qubit)
    for qubit in annhiliation_qubit_string:
        if qubit not in creation_qubit_string:
            nonoverlapping_qubits.append(qubit)
    return nonoverlapping_qubits

def count_pauli_weight_for_interaction(creation_mode_coordinates, annhiliation_mode_coordinates, lattice: Lattice, fenwick_trees):
    creation_update_qubits, creation_target_qubit, creation_parity_qubits = gates_needed_to_create_annhilate_fermionic_mode(lattice, creation_mode_coordinates, fenwick_trees)
    annhiliation_update_qubits, annhiliation_target_qubit, annhiliation_parity_qubits = gates_needed_to_create_annhilate_fermionic_mode(lattice, annhiliation_mode_coordinates, fenwick_trees)
    non_overlapping_update_qubits = find_nonoverlapping_qubits(creation_update_qubits,annhiliation_update_qubits)
    for qubit in non_overlapping_update_qubits:
        print("Update qubit: ", qubit.root_enumeration, ", " ,qubit.tree_enumeration)
    non_overlapping_target_qubits = [annhiliation_target_qubit, creation_target_qubit] if annhiliation_target_qubit != creation_target_qubit else []
    non_overlapping_parity_qubits = find_nonoverlapping_qubits(creation_parity_qubits, annhiliation_parity_qubits)
    for qubit in non_overlapping_parity_qubits:
        print("Parity qubit: ", qubit.root_enumeration, ", ", qubit.tree_enumeration)
    non_overlapping_qubits = non_overlapping_parity_qubits+ non_overlapping_target_qubits  + non_overlapping_update_qubits
    return len(set(non_overlapping_qubits))

def construct_fenwick_tree(parent, root_enumeration, L, R) -> List[Qubit]:
    if(L != R):
        qubit = Qubit(parent, root_enumeration)
        if parent != None:
            parent.add_child(qubit)
        fenwick_tree = [qubit]
        fenwick_tree.extend(construct_fenwick_tree(qubit, root_enumeration, math.floor((R+L)/2) +1, R))
        fenwick_tree.extend(construct_fenwick_tree(qubit, root_enumeration, L, math.floor((R+L)/2)))
        return fenwick_tree
    else:
        return []

def enumerate_fenwick_tree(fenwick_tree):
    for i, qubit in enumerate(fenwick_tree):
        qubit.set_tree_enumeration(len(fenwick_tree) - i - 1)
    return fenwick_tree

def construct_fenwick_trees(size_of_tree, number_of_trees) -> List[List[Qubit]]:
    fenwick_trees = []
    for i in range(number_of_trees):
        fenwick_tree = construct_fenwick_tree(parent=None,root_enumeration=i, L=0, R=size_of_tree)
        fenwick_tree = enumerate_fenwick_tree(fenwick_tree)
        fenwick_trees.append(fenwick_tree)
    return fenwick_trees

def fenwick_trees_to_qubit_list(fenwick_trees: List[List[Qubit]]):
    qubits = []
    for fenwick_tree in fenwick_trees:
        for qubit in fenwick_tree:
            qubits.append(qubit)
    return qubits

def map_fenwick_trees_to_lattice(fenwick_trees: List[List[Qubit]],width) -> Lattice:
    num_of_trees = len(fenwick_trees)
    size_of_trees = len(fenwick_trees[0])
    max_size_of_lattice = num_of_trees*size_of_trees
    number_of_rows = max_size_of_lattice // width
    qubits = fenwick_trees_to_qubit_list(fenwick_trees)
    if (size_of_trees == 4 and width % 4 == 0):
        enumeration_x_coordinates, enumeration_y_coordinates = map_fenwick_trees_to_lattice_four(fenwick_trees, width)
        lattice = Lattice(number_of_rows, width,qubits, enumeration_x_coordinates, enumeration_y_coordinates)
        return lattice

def map_fenwick_trees_to_lattice_four(fenwick_trees, width):
    enumeration_x_coordinates = []
    enumeration_y_coordinates = []
    for fenwick_tree in fenwick_trees:
         enumeration_x_coordinates.append(fenwick_tree[0].root_enumeration*2 % width)
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width) *2)
         enumeration_x_coordinates.append(fenwick_tree[0].root_enumeration*2 % width)
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width)*2+1)
         enumeration_x_coordinates.append(fenwick_tree[0].root_enumeration*2 % width +1)
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width)*2)
         enumeration_x_coordinates.append(fenwick_tree[0].root_enumeration*2 % width+1)
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width) *2+1)
    return enumeration_x_coordinates, enumeration_y_coordinates
size_of_fenwick_trees = 4
N= 20
number_of_nodes = N**2
width_of_lattice = N
number_of_trees = math.ceil(number_of_nodes/ size_of_fenwick_trees)
fenwick_trees = construct_fenwick_trees(size_of_fenwick_trees,number_of_trees)

lattice = map_fenwick_trees_to_lattice(fenwick_trees,width_of_lattice)
total_pauli_weight = 0
for coordinate in lattice.coordinates:
        neighbours_coordinates = lattice.get_neighbours(coordinate[0],coordinate[1])
        for neighbour_coordinate in neighbours_coordinates:
            print(coordinate," - ",neighbour_coordinate)
            pauli_weight = count_pauli_weight_for_interaction(coordinate,neighbour_coordinate,lattice,fenwick_trees)
            print(pauli_weight)
            total_pauli_weight += pauli_weight
number_of_interactions = (lattice.width - 1) * lattice.height + (lattice.height - 1)* lattice.width
print(number_of_interactions)
print(total_pauli_weight)
print((N**3 - N)/number_of_interactions + 1)
print(total_pauli_weight/(2*number_of_interactions))
#lattice = Lattice(height=5, width=5, qubits=qubits,enumeration_x_coordinates=enumeration_x_coordinates, enumeration_y_coordinates=enumeration_y_coordinates)