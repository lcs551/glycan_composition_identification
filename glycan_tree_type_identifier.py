# Identify whether a glycan tree is high mannose, complex or hybrid.
# python glycan_tree_type_identified.py --wurcs "WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/"

import re
import argparse
from typing import List

database = {
 "a2122h-1a_1-5":"GLC" ,
 "a2122h-1b_1-5":"BGC" , 
 "a1122h-1a_1-5":"MAN" , 
 "a1122h-1b_1-5":"BMA" ,
 "a2112h-1a_1-5":"GLA" , 
 "a2112h-1b_1-5":"GAL" , 
 "a2112m-1a_1-5":"FCA" , 
 "a2112m-1b_1-5":"FCB" ,
 "a1221m-1a_1-5":"FUC" , 
 "a1221m-1b_1-5":"FUL" , 
 "a212h-1a_1-5" :"XYS" , 
 "a212h-1b_1-5" :"XYP" , 
 "a21d2h-1b_1-5":"Z9D" ,
 "a2122h-1b_1-5_2*N":"GCS",
 "a2122h-1a_1-5_2*N":"PA1",
 "a2122h-1b_1-5_2*NCC/3=O": "NAG",
 "a2122h-1a_1-5_2*NCC/3=O": "NDG",
 "a2112h-1b_1-5_2*NCC/3=O": "NGA",
 "a2112h-1a_1-5_2*NCC/3=O": "A2G",
 "a1122h-1a_1-5_2*NCC/3=O": "BM3",
 "a1122h-1b_1-5_2*NCC/3=O": "BM7"
 }

def get_unique_sugars(WURCS: str):
    """
    Find unique sugars in the given WURCS string.

    :param WURCS: The WURCS string to search for unique sugars.
    :return: A list of sugar names corresponding to the unique sugars found in the WURCS string.
    """
    if 'ERROR' in WURCS:
        return
    sugar_regex ="\[\S*\]"
    sugars = re.findall(sugar_regex, WURCS)
    x = sugars[0]
    x = x.lstrip("[").rstrip("]")
    y = x.split("][")

    sugar_names = []

    # Search database
    for string in y:
        if string in database:
            sugar_names.append(database[string])
        else:
            sugar_names.append(None)

    return sugar_names

def get_sugar_order(WURCS: str):
    """
    :param WURCS: The WURCS code of the sugar molecule.
    :return: The list of sugar order. Each element represents the order of a sugar.
    """
    order_regex = "/[0-9-]*/"
    order = re.findall(order_regex, WURCS)[0]
    return order.replace("/","").split("-")

def get_linkages(WURCS: str):
    """
    :param WURCS: The WURCS string from which linkages need to be extracted.
    :return: A list of linkages found in the WURCS string.
    """
    linkage_regex = "[a-z]{1}[0-9]{1}-{1}[a-z]{1}[0-9]"
    linkages = re.findall(linkage_regex, WURCS)
    return linkages

def organise_linkages(linkages: List[str]):
    """
    Organises linkages into branches based on a specific branch point.

    :param linkages: List of linkages in the format 'donor-acceptor'.
    :return: List of branches, where each branch is a list of nodes.
    """
    links = {}

    for linkage in linkages:
        donor, acceptor = linkage.split("-")
        links.setdefault(donor[0], []).append(acceptor[0])

    ### Added count so the branchpoint is set as the first node ###
    ### rather than last node, which lead to branches being wrong ###    
    for k, v in links.items():
        if len(v) > 1:
            branchpoint = k
            break
    else:
        # This part will be executed if the loop completes without a break
        return

    def search(node, curr): 
        if node not in links: 
            curr.append(node)
            return
        else:
            curr.append(node)
        
        for child in links[node]:
            search(child, curr)
            
        return curr

    branches = []
    for node in links[branchpoint]:
        output = []
        search(node, output)
        branches.append(output)
    
    return branches

def branches_to_sugars(branches: List, sugar_alphabet_map: dict):
    """
    Convert branches to sugars based on the given sugar alphabet map.

    :param branches: A list of branches, where each branch is a list of alphabets.
    :param sugar_alphabet_map: A dictionary mapping alphabets to corresponding sugars.
    :return: A list of converted sugar branches.

    Example:
    >>> branches = [['A', 'B', 'C'], ['D', 'E', 'F']]
    >>> sugar_alphabet_map = {'A': 'SugarA', 'B': 'SugarB', 'C': 'SugarC', 'D': 'SugarD', 'E': 'SugarE', 'F': 'SugarF'}
    >>> branches_to_sugars(branches, sugar_alphabet_map)
    [['SugarA', 'SugarB', 'SugarC'], ['SugarD', 'SugarE', 'SugarF']]
    """
    sugar_branches = []
    for branch in branches:
        tmp = []
        for alphabet in branch:
            tmp.append(sugar_alphabet_map[alphabet])
        sugar_branches.append(tmp)    
    return sugar_branches

def check_type(WURCS: str):
    """
    The `check_type` method is used to determine the type of a glycan based on its WURCS string representation.

    :param WURCS: The WURCS string representation of the glycan.
    :return: The type of the glycan, which can be "High Mannose", "Hybrid", or "Complex".

    """
    sugars = get_unique_sugars(WURCS=WURCS)
    if sugars == None:
        return "Error producing WURCS string"
    order = get_sugar_order(WURCS=WURCS)
    sugar_list = [sugars[int(num) - 1] for num in order] # Correspond sugar names to their order
    suitable_glycan = 0

    # Check if there is a suitable glycan core: 
    # Must have MAN/BMA residue to be long enough to be considered (excludes glycan chains of just NAG or NAG, NAG)
    # Also excludes e.g. 6DTU where 'glycan' cain is GLC,GLC,GLC,GLC
    for sugar in sugar_list:
        if sugar == 'MAN' or sugar == 'BMA':
            suitable_glycan +=1            
    
    if suitable_glycan == 0:
        return "Unsuitable core glycan"

    linkages = get_linkages(WURCS=WURCS)
    branches = organise_linkages(linkages=linkages)

    sugar_map = {}
    alphabet_map = {}
    alphabet = "abcdefghijklmnopqrstuvxyz"

    for index, pos in enumerate(order): 
        pos = int(pos)
        sugar_map[pos] = sugars[pos-1]
        alphabet_map[alphabet[index]] = pos

    sugar_alphabet_map = {}
    for k, v in alphabet_map.items(): 
        sugar_alphabet_map[k] = sugar_map[v]

    ### Identify whether the glycan tree is high mannose or not, by seeing if any ###
    ### MAN/BMA residues are found in the list after the first MAN residue is found ###
    found_man = False

    for value in sugar_list:
        if found_man:
            if value != 'MAN' and value != 'BMA':
                found_man = False
                break
        elif value == 'MAN' or value == 'BMA':
            found_man = True

    if found_man:
        return "High Mannose"
   
    
    branches = branches_to_sugars(branches=branches, sugar_alphabet_map=sugar_alphabet_map)    

    # Check how many of the branches are mannose only, if any of them are, then it is a hybrid otherwise it is complex
        
    mannose_only_branches = 0
    for branch in branches:
        if list(set(branch)) == ["MAN"]:
            mannose_only_branches+=1
    
    if mannose_only_branches > 0:
        return "Hybrid"

    return "Complex"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="glycan_composition_identification",
        description="""Identify whether your glycan tree is high mannose, complex or hybrid."""
    )
    parser.add_argument(
        "-w ",
        "--wurcs"
    )

    args = parser.parse_args()
    user_wurcs = args.wurcs

    hm = "WURCS=2.0/3,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-e1_e3-f1_f2-g1"
    complex = "WURCS=2.0/3,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-1-3-1/a4-b1_b4-c1_c3-d1_c6-f1_d2-e1_f2-g1"
    hyb= "WURCS=2.0/4,8,7/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5]/1-2-3-1-4-3-3-3/a4-b1_b3-c1_b6-f1_c2-d1_d4-e1_f3-g1_f6-h1"
    hm_branched = "WURCS=2.0/3,11,10/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3-3-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-g1_d2-e1_e2-f1_g3-h1_g6-j1_h2-i1_j2-k1"
    hm_linear = "WURCS=2.0/3,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3/a4-b1_b4-c1_c3-d1_d2-e1_e2-f1"

    unsuitable_one_NAG = "WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/" 
    unsuitable_one_GLC = "WURCS=2.0/1,4,3/[a2122h-1a_1-5]/1-1-1-1/a4-b1_b4-c1_c4-d1" # GLC, GLC, GLC, GLC

    branch_issue_1 = "WURCS=2.0/4,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5][a1122h-1b_1-5][a1122h-1a_1-5]/1-2-1-3-4-4-2/a3-b1_a4-c1_a6-g1_c4-d1_d3-e1_d6-f1"
    branch_issue_2 = "WURCS=2.0/3,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-1-3/a4-b1_b4-c1_c3-d1_c6-f1_d2-e1"

    result = check_type("WURCS=2.0/3,3,2/[a1122h-1a_1-5_1*OC][a1122h-1a_1-5][a2122h-1b_1-5_2*NCC/3=O]/1-2-3/a3-b1_b2-c1"
)
    print(result)