from openbabel import openbabel as ob
import numpy as np
from numpy import linalg
import math

class OBProtein():

    @staticmethod
    def Hydrogenate(carbon, mainAtom, mol):

        """
        Adds 3 hydrogens to a molecule that will be bonded to a carbon

        carbon: The carbon atom where the 3 hydrogens will be added to
        mainAtom: The other atom the carbon is attached to
        mol: The molecule that will be hydrogenated

        returns: The hydrogenated molecule
        """

        chbl = 1.09

        bond_length = carbon.GetDistance(mainAtom)

        Ux = mainAtom.GetX() - carbon.GetX()
        Uy = mainAtom.GetY() - carbon.GetY()
        Uz = mainAtom.GetZ() - carbon.GetZ()

        U = np.array([Ux, Uy, Uz])

        V = np.array([Uy*Uz, -2*Ux*Uz, Ux*Uy])

        U = U/(np.linalg.norm(U))

        V = V/(np.linalg.norm(V))

        ang1 = math.radians(70.5)

        H1 = -U*math.cos(ang1) + V*math.sin(ang1)

        h1 = ob.OBAtom()
        h1.SetAtomicNum(1)
        h1.SetVector(carbon.GetX() + H1[0]*chbl, carbon.GetY() + H1[1]*chbl, carbon.GetZ() + H1[2]*chbl)

        mol.AddAtom(h1)

        Z = -U - H1
        Z = Z/(np.linalg.norm(Z))

        Y = np.cross(U, H1)
        Y = Y/(np.linalg.norm(Y))

        X = np.cross(Y, Z)
        X = X/(np.linalg.norm(X))

        ang2 = math.radians(54.7)

        H2 = Z*math.cos(ang2) + Y*math.sin(ang2)

        h2 = ob.OBAtom()
        h2.SetAtomicNum(1)

        h2.SetVector(carbon.GetX() + H2[0]*chbl, carbon.GetY() + H2[1]*chbl, carbon.GetZ() + H2[2]*chbl)

        if mol.AddAtom(h2):
            print("Success adding H2")

        H3 = -(H1 + H2 + U)
        H3 = H3/(np.linalg.norm(H3))

        h3 = ob.OBAtom()
        h3.SetAtomicNum(1)

        h3.SetVector(carbon.GetX() + H3[0]*chbl, carbon.GetY() + H3[1]*chbl, carbon.GetZ() + H3[2]*chbl)

        if mol.AddAtom(h3):
            print("Success adding H3")

        return mol

    @staticmethod
    def IsInMolecule(atom, molecule):
        for mol_atom in ob.OBMolAtomIter(molecule):
            if OBAminoAcid.Equals(mol_atom, atom):
                return True
        return False

    @staticmethod
    def Equals(atom1, atom2):
        return (atom1.GetX() == atom2.GetX()) and (atom1.GetY() == atom2.GetY()) and (atom1.GetZ() == atom2.GetZ()) and (atom1.GetAtomicNum() == atom2.GetAtomicNum())

    def __init__(self, mol):
        self.mol = mol
        self.SetResiduesPerceived()

    def SetMolecule(self, mol):
        self.mol = mol
        self.SetResiduesPerceived()

    def SetResiduesPerceived(self):
        self.residues = []

        count = 0
        for atom1 in ob.OBMolAtomIter(self.mol):
            if atom1.GetAtomicNum() == 7:
                for atom2 in ob.OBAtomAtomIter(atom1):
                    if atom2.GetAtomicNum() == 6:
                        for atom3 in ob.OBAtomAtomIter(atom2):
                            if atom3.GetAtomicNum() == 6:
                                for atom4 in ob.OBAtomAtomIter(atom3):
                                    #print(atom4.GetAtomicNum(), atom3.GetBond(atom4).GetBondOrder())
                                    if atom4.GetAtomicNum() == 8 and atom3.GetBond(atom4).GetBondOrder() == 2:
                                        count += 1
                                        self.residues.append(OBAminoAcid(self, count, atom1, atom3))

        #print(count)

    def CapAllResidues(self):

        for i in range(len(self.residues)):
            res = self.residues[i]
            resMol = res.resMol
            nTerm = res.NTerminus
            cTerm = res.CTerminus

            nTermIsCapped = False
            cTermIsCapped = False

            for atom1 in ob.OBAtomAtomIter(nTerm):
                if OBAminoAcid.Equals(nTerm, atom1):
                    continue
                elif atom1.GetAtomicNum() == 6:
                    for atom2 in ob.OBAtomAtomIter(atom1):
                        if OBAminoAcid.Equals(nTerm, atom2):
                            continue
                        elif atom2.GetAtomicNum() == 8 and atom1.GetBond(atom2).GetBondOrder() == 2:
                            self.residues[i].resMol.AddAtom(atom1)
                            self.residues[i].resMol.AddAtom(atom2)
                            self.residues[i].NTermCap.AddAtom(atom1)
                            self.residues[i].NTermCap.AddAtom(atom2)
                            a1 = atom1
                            nTermIsCapped = True

            if nTermIsCapped:
                for atom3 in ob.OBAtomAtomIter(a1):
                    if atom3.GetAtomicNum() == 6:
                        self.residues[i].resMol.AddAtom(atom3)
                        self.residues[i].NTermCap.AddAtom(atom3)
                        self.residues[i].resMol = OBProtein.Hydrogenate(atom3, atom1, self.residues[i].resMol)
                        self.residues[i].NTermCap = OBProtein.Hydrogenate(atom3, atom1, self.residues[i].NTermCap)
                        break

            if not nTermIsCapped:
                for atom in ob.OBAtomAtomIter(nTerm):
                    if not IsInMolecule(atom, self.residues[i].resMol):
                        self.residues[i].resMol.AddAtom(atom)

            for atom1 in ob.OBAtomAtomIter(cTerm):
                if atom1.GetAtomicNum() == 7:
                    for atom2 in ob.OBAtomAtomIter(atom1):
                        if not OBProtein.Equals(cTerm, atom2):
                            a1 = atom1
                            self.residues[i].resMol.AddAtom(atom2)
                            self.residues[i].CTermCap.AddAtom(atom2)
                            if atom2.GetAtomicNum() == 6:
                                self.residues[i].resMol = OBProtein.Hydrogenate(atom2, atom1, self.residues[i].resMol)
                                self.residues[i].CTermCap = OBProtein.Hydrogenate(atom2, atom1, self.residues[i].CTermCap)
                                cTermIsCapped = True

            if cTermIsCapped:
                self.residues[i].resMol.AddAtom(a1)
                self.residues[i].CTermCap.AddAtom(a1)

            if not cTermIsCapped:
                for atom in ob.OBAtomAtomIter(cTerm):
                    if atom.GetAtomicNum() == 8 and atom.GetBond(cTerm).GetBondOrder() == 1:
                        self.residues[i].resMol.AddAtom(atom)

class OBAminoAcid():

    @staticmethod
    def Equals(atom1, atom2):
        return (atom1.GetX() == atom2.GetX()) and (atom1.GetY() == atom2.GetY()) and (atom1.GetZ() == atom2.GetZ()) and (atom1.GetAtomicNum() == atom2.GetAtomicNum())

    @staticmethod
    def IsInMolecule(atom, molecule):
        for mol_atom in ob.OBMolAtomIter(molecule):
            if OBAminoAcid.Equals(mol_atom, atom):
                return True
        return False

    @staticmethod
    def SideChainGrowerRecursive(atom, mol, cmol):
        for n_atom in ob.OBAtomAtomIter(atom):
            print(atom.GetAtomicNum(), n_atom.GetAtomicNum(), OBAminoAcid.IsInMolecule(n_atom, mol))
            if not OBAminoAcid.IsInMolecule(n_atom, mol):
                mol.AddAtom(n_atom)
                cmol.AddAtom(n_atom)
                OBAminoAcid.SideChainGrower(n_atom, mol, cmol)

    @staticmethod
    def SideChainGrowerHelper(atom, mol, cmol):
        count = 0
        for n_atom in ob.OBAtomAtomIter(atom):
            print(atom.GetAtomicNum(), n_atom.GetAtomicNum(), OBAminoAcid.IsInMolecule(n_atom, mol))
            if not OBAminoAcid.IsInMolecule(n_atom, mol) and not (atom.GetAtomicNum() == 16 and n_atom.GetAtomicNum() == 16):
                mol.AddAtom(n_atom)
                cmol.AddAtom(n_atom)
                count += 1

        return count

    @staticmethod
    def SideChainGrowerIterative(protein, mol, cmol):
        count = 1
        while count > 0:
            count = 0
            for mol_atom in ob.OBMolAtomIter(protein):
                if OBAminoAcid.IsInMolecule(mol_atom, cmol):
                    count += OBAminoAcid.SideChainGrowerHelper(mol_atom, mol, cmol)

    @staticmethod
    def BondOrderSum(atom):
        bosum = 0
        for n_atom in ob.OBAtomAtomIter(atom):
            bosum += atom.GetBond(n_atom).GetBondOrder()

        return bosum

    @staticmethod
    def SulfurHydrogenator(sulfur, mainAtom, mol):

        sulfur_coords = np.array([sulfur.GetX(), sulfur.GetY(), sulfur.GetZ()], dtype=float)
        main_atom_coords = np.array([mainAtom.GetX(), mainAtom.GetY(), mainAtom.GetZ()], dtype=float)

        S = sulfur_coords
        C = main_atom_coords

        V = C - S

        print(S, C, V)

        V = V/np.linalg.norm(V)

        R = np.array([V[1]*V[2], -2*V[0]*V[2], V[0]*V[1]], dtype=float)

        R = R/np.linalg.norm(R)

        H = R*math.cos(math.radians(2.1)) - V*math.sin(math.radians(2.1))

        print(V, R, H)

        H = S + (1.34 * H)

        print(S, H)

        hydrogen = ob.OBAtom()
        hydrogen.SetAtomicNum(1)

        print(H[0], H[1], H[2])

        hydrogen.SetVector(H[0], H[1], H[2])

        if mol.AddAtom(hydrogen):
            print("Hydrogen Added Successfully to Sulfur")

        return mol

    def __init__(self, protein, residueNumber, NTerminus, CTerminus):
        self.protein = protein
        self.residueNumber = residueNumber
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus
        self.resMol = ob.OBMol()
        self.sideChain = ob.OBMol()
        self.NTermCap = ob.OBMol()
        self.CTermCap = ob.OBMol()

        self.resMol.AddAtom(NTerminus)
        self.resMol.AddAtom(CTerminus)

        for atomN in ob.OBAtomAtomIter(NTerminus):
            for atomC in ob.OBAtomAtomIter(CTerminus):
                if OBAminoAcid.Equals(atomN, atomC):
                    self.AlphaCarbon = atomC
                    self.resMol.AddAtom(atomC)

        for atom in ob.OBAtomAtomIter(NTerminus):
            if atom.GetAtomicNum() == 1:
                self.resMol.AddAtom(atom)

        for atom in ob.OBAtomAtomIter(CTerminus):
            if atom.GetAtomicNum() == 8 and atom.GetBond(CTerminus).GetBondOrder() == 2:
                self.resMol.AddAtom(atom)

        for atom in ob.OBAtomAtomIter(self.AlphaCarbon):
            if atom.GetAtomicNum() == 1:
                self.resMol.AddAtom(atom)
                self.AlphaHydrogen = atom
                break

        for atom in ob.OBAtomAtomIter(self.AlphaCarbon):
            if not (OBAminoAcid.Equals(atom, self.NTerminus) or OBAminoAcid.Equals(atom, self.CTerminus) or OBAminoAcid.Equals(atom, self.AlphaHydrogen)):
                self.resMol.AddAtom(atom)
                self.sideChain.AddAtom(atom)
                print(self.residueNumber)

        OBAminoAcid.SideChainGrowerIterative(self.protein.mol, self.resMol, self.sideChain)

        hasDisulfide = False

        for atom in ob.OBMolAtomIter(self.protein.mol):
            print(F'{atom.GetAtomicNum()},BOSum: {OBAminoAcid.BondOrderSum(atom)}')
            if atom.GetAtomicNum() == 16 and OBAminoAcid.IsInMolecule(atom, self.resMol):
                for n_atom in ob.OBAtomAtomIter(atom):
                    if n_atom.GetAtomicNum() == 16:
                        hasDisulfide = True

        if hasDisulfide:
            for atom in ob.OBMolAtomIter(self.protein.mol):
                if atom.GetAtomicNum() == 16 and OBAminoAcid.IsInMolecule(atom, self.resMol):
                    for n_atom in ob.OBAtomAtomIter(atom):
                        if n_atom.GetAtomicNum() == 6:
                            self.resMol = OBAminoAcid.SulfurHydrogenator(atom, n_atom, self.resMol)
                            self.sideChain = OBAminoAcid.SulfurHydrogenator(atom, n_atom, self.sideChain)

        #Set the Residue Name of the Amino Acid

        countH = 0
        countC = 0
        countN = 0
        countO = 0
        countS = 0

        for atom in ob.OBMolAtomIter(self.sideChain):
            if atom.GetAtomicNum() == 1:
                countH += 1
            elif atom.GetAtomicNum() == 6:
                countC += 1
            elif atom.GetAtomicNum() == 7:
                countN += 1
            elif atom.GetAtomicNum() == 8:
                countO += 1
            elif atom.GetAtomicNum() == 16:
                countS += 1

        countTuple = (countH, countC, countN, countO, countS)

        if countTuple == (1, 0, 0, 0, 0):
            self.name = "GLY"
        elif countTuple == (3, 1, 0, 0, 0):
            self.name = "ALA"
        elif countTuple == (7, 3, 0, 0, 0):
            self.name = "VAL"
        elif countTuple == (9, 4, 0, 0, 0):
            carbonCount = 0
            for atom in ob.OBAtomAtomIter(self.AlphaCarbon):
                if OBAminoAcid.IsInMolecule(atom, self.SideChain):
                    for n_atom in ob.OBAtomAtomIter(atom):
                        if n_atom.GetAtomicNum() == 6:
                            carbonCount += 1

            if carbonCount == 2:
                self.name = "LEU"
            elif carbonCount == 3:
                self.name = "ILE"

        elif countTuple == (7, 3, 0, 0, 1):
            self.name = "MET"
        elif countTuple == (7, 7, 0, 0, 0):
            self.name = "PHE"
        elif countTuple == (8, 9, 1, 0, 0):
            self.name = "TRP"
        elif countTuple == (6, 3, 0, 0, 0):
            self.name = "PRO"
        elif countTuple == (3, 1, 0, 1, 0):
            self.name = "SER"
        elif countTuple == (5, 2, 0, 1, 0):
            self.name = "THR"
        elif countTuple == (7, 7, 0, 1, 0):
            self.name = "TYR"
        elif countTuple == (3, 1, 0, 0, 1):
            self.name = "CYS"
        elif countTuple == (2, 2, 0, 2, 0):
            self.name = "ASP"
        elif countTuple == (11, 4, 1, 0, 0):
            self.name = "LYS"
        elif countTuple == (4, 2, 1, 1, 0):
            self.name = "ASN"
        elif countTuple == (6, 3, 1, 1, 0):
            self.name = "GLN"
        elif countTuple == (4, 3, 0, 2, 0):
            self.name = "GLU"
        elif countTuple == (11, 4, 3, 0, 0):
            self.name = "ARG"
        elif countTuple == (5, 4, 2, 0, 0) or countTuple == (6, 4, 2, 0, 0):
            self.name = "HIS"

class OBLigand():

    @staticmethod
    def Equals(atom1, atom2):
        return (atom1.GetX() == atom2.GetX()) and (atom1.GetY() == atom2.GetY()) and (atom1.GetZ() == atom2.GetZ()) and (atom1.GetAtomicNum() == atom2.GetAtomicNum())

    @staticmethod
    def IsInMolecule(atom, molecule):
        for mol_atom in ob.OBMolAtomIter(molecule):
            if OBAminoAcid.Equals(mol_atom, atom):
                return True
        return False

    @staticmethod
    def IsConnectedTo(atom, mol, big_mol):

        if OBLigand.IsInMolecule(atom, mol):
            return False

        for mol_atom in ob.OBMolAtomIter(big_mol):
            if OBLigand.IsInMolecule(mol_atom, mol):
                for n_atom in ob.OBAtomAtomIter(mol_atom):
                    if OBLigand.Equals(n_atom, atom):
                        return True

        return False

    @staticmethod
    def MolGrowerHelper(atom, mol):
        count = 0
        for n_atom in ob.OBAtomAtomIter(atom):
            #print(atom.GetAtomicNum(), n_atom.GetAtomicNum(), OBLigand.IsInMolecule(n_atom, mol))
            if not OBLigand.IsInMolecule(n_atom, mol):
                mol.AddAtom(n_atom)
                count += 1

        return count

    @staticmethod
    def MolGrowerIterative(mol, big_mol):
        atomArr = []
        count = 1
        while count > 0:
            count = 0
            for mol_atom in ob.OBMolAtomIter(big_mol):
                if OBLigand.IsInMolecule(mol_atom, mol):
                    count += OBLigand.MolGrowerHelper(mol_atom, mol)

    @staticmethod
    def SplitUp(mol):
        molArr = []
        addNew = True

        for atom in ob.OBMolAtomIter(mol):
            for i in range(len(molArr)):
                if OBLigand.IsInMolecule(atom, molArr[i]):
                    addNew = False

            if addNew:
                molArr.append(ob.OBMol())
                molArr[-1].AddAtom(atom)
                OBLigand.MolGrowerIterative(molArr[-1], mol)

            addNew = True

        return molArr

    @staticmethod
    def IsProtein(mol, big_mol):

        for atom1 in ob.OBMolAtomIter(big_mol):
            if OBLigand.IsInMolecule(atom1, mol):
                if atom1.GetAtomicNum() == 7:
                    for atom2 in ob.OBAtomAtomIter(atom1):
                        if atom2.GetAtomicNum() == 6:
                            for atom3 in ob.OBAtomAtomIter(atom2):
                                if atom3.GetAtomicNum() == 6:
                                    for atom4 in ob.OBAtomAtomIter(atom3):
                                        if atom4.GetAtomicNum() == 8 and atom3.GetBond(atom4).GetBondOrder() == 2:
                                            return True
        return False


    def __init__(self, whole_mol):
        self.whole_mol = whole_mol
        self.ligand_mols = []

        molArr = OBLigand.SplitUp(self.whole_mol)

        for i in range(len(molArr)):
            if not OBLigand.IsProtein(molArr[i], self.whole_mol):
                self.ligand_mols.append(molArr[i])
