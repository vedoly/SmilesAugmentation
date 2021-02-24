import random
from rdkit import *
import numpy as np

count = 0


def randomSmiles(m1, n):
    global count
    count += 1
    s = set()
    if type(m1) == str:
        m1 = Chem.MolFromSmiles(m1)
    for i in range(n):
        m1.SetProp("_canonicalRankingNumbers", "True")
        idxs = list(range(0, m1.GetNumAtoms()))
        random.shuffle(idxs)
        for i, v in enumerate(idxs):
            m1.GetAtomWithIdx(i).SetProp("_canonicalRankingNumber", str(v))
        s.add(Chem.MolToSmiles(m1))
    print(count)
    return s


def randomSmilesExpand(m1, start_n=2000, multiplier=10):
    first = set()
    second = set()
    for i in range(100):
        first = first | randomSmiles(m1, start_n * (multiplier ** i))
        second = second | randomSmiles(m1, start_n * (multiplier ** i))
        # print(len(first), len(second))
        if len(first) == len(second):
            return first
        # print(i + 1)


def randomize_smiles(smiles):
    """Perform a randomization of a SMILES string
        must be RDKit sanitizable"""
    m = Chem.MolFromSmiles(smiles)
    ans = list(range(m.GetNumAtoms()))
    np.random.shuffle(ans)
    nm = Chem.RenumberAtoms(m, ans)
    return nm


# m1 = Chem.MolFromSmiles("CC(C)(C)OC(=O)OC(=O)OC(C)(C)C")
# s = randomSmiles(m1,30)
