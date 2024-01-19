import re


WURCS = "WURCS=2.0/3,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-e1_e3-f1_f2-g1"


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
    order_regex = "/[0-9-]*/"
    order = re.findall(order_regex, WURCS)[0]
    return order.replace("/","").split("-")

def get_linkages(WURCS: str):
    linkage_regex = "[[a-z]{1}[0-9]{1}-{1}[a-z]{1}[0-9]"
    linkages = re.findall(linkage_regex, WURCS)
    return linkages

def main(): 
    sugars = get_unique_sugars(WURCS=WURCS)
    linkages = get_linkages(WURCS=WURCS)
    order = get_sugar_order(WURCS=WURCS)

    sugar_map = {}
    alphabet_map = {}
    alphabet = "abcdefghijklmnopqrstuvxyz"

    for index, pos in enumerate(order): 
        pos = int(pos)
        sugar_map[pos] = sugars[pos-1]
        alphabet_map[alphabet[index]] = pos

    print(sugar_map)
    print(alphabet_map)

    ### Correspond sugar names to their order ###
    sugar_list = [sugars[int(num) - 1] for num in order]
    # print(sugar_list)

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
        result = "High Mannose"
    else:
        result = "Not High Mannose"

    # print(result)
    return result
    
    # if all(string in ['NAG', 'MAN'] for string in result): # can't do this because some hybrid glycans have NAG at end of chain
    #     print("High mannose")
    # else:
    #     print('Not high mannose')

if __name__ == "__main__":
    main()