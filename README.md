
The Glycan Tree Type Identifier identifies the tree type of an N-glycan using the WURCS code.
There are four possible outputs:
  - High mannose
  - Complex
  - Hybrid
  - Unsuitable core glycan

Unsuitable core glycan represents a glycan chain which is too short to designate as high mannose, complex or hybrid.

The most up to date scripts are in the /scripts folder.
To run the glycan composition identifier: python glycan_tree_type_identifier.py --wurcs "WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/"
