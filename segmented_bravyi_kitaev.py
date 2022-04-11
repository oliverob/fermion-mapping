
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
    
    def set_tree_enumeration(self, tree_enumeration):
        self.tree_enumeration = tree_enumeration

    def get_update_qubits(self):
        if self.parent == None:
            return None
        return [self.parent, self.get_parent().get_update_qubits()]


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
            left_neighbour = self.lattice_array[x-1, y]
            neighbours.append(left_neighbour)
        if(x < self.width):
            right_neigbour = self.lattice_array[x+1, y]
            neighbours.append(right_neigbour)
        if(y > 0):
            above_neighbour = self.lattice_array[x, y-1]
            neighbours.append(above_neighbour)
        if(y < self.height):
            below_neighbour = self.lattice_array[x, y+1]
            neighbours.append(below_neighbour)
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
    update_qubits = target_qubit.get_update_qubits()
    parity_qubits = fenwick_trees[target_qubit.root_enumeration][:target_qubit.tree_enumeration]
    parity_qubits.extend([fenwick_trees[i][0] for i in range(0,target_qubit.tree_enumeration)])
    return [target_qubit, update_qubits, parity_qubits]

def count_pauli_weight_for_interaction(creation_mode_coordinates, annhiliation_mode_coordinates, lattice: Lattice, fenwick_trees):
    creation_qubits = gates_needed_to_create_annhilate_fermionic_mode(lattice, creation_mode_coordinates, fenwick_trees)
    annhiliation_qubits = gates_needed_to_create_annhilate_fermionic_mode(lattice, annhiliation_mode_coordinates, fenwick_trees)
    non_overlapping_qubits = [qubit for qubit in creation_qubits if qubit not in annhiliation_qubits]
    return len(non_overlapping_qubits)

def construct_fenwick_tree(parent, root_enumeration, L, R) -> List[Qubit]:
    if(L != R):
        qubit = Qubit(parent, root_enumeration)
        fenwick_tree = [qubit]
        fenwick_tree.extend(construct_fenwick_tree(qubit, root_enumeration, L, math.floor((R+L)/2)))
        fenwick_tree.extend(construct_fenwick_tree(qubit, root_enumeration, math.floor((R+L)/2) +1, R))
        return fenwick_tree
    else:
        return []

def enumerate_fenwick_tree(fenwick_tree):
    for i, qubit in enumerate(fenwick_tree):
        qubit.set_tree_enumeration(i)
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
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width)*2+1)
         enumeration_x_coordinates.append(fenwick_tree[0].root_enumeration*2 % width+1)
         enumeration_y_coordinates.append((fenwick_tree[0].root_enumeration*2 // width) *2)
    return enumeration_x_coordinates, enumeration_y_coordinates

fenwick_trees = construct_fenwick_trees(4,8)
print(fenwick_trees)
print(fenwick_trees[0][2].parent)
lattice = map_fenwick_trees_to_lattice(fenwick_trees,4)
print(count_pauli_weight_for_interaction([0,1],[0,0],lattice,fenwick_trees))
#lattice = Lattice(height=5, width=5, qubits=qubits,enumeration_x_coordinates=enumeration_x_coordinates, enumeration_y_coordinates=enumeration_y_coordinates)