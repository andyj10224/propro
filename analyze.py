#Andy Jiang, Sherrill Group, Georgia Institute of Technology

from openbabel import openbabel as ob
import numpy as np
from numpy import linalg
import math
import sys
import subprocess
import os

def GetCharge(mol):
    """
    Gets the charge of an OBMol object (Only works for protein residues)

    mol: An OBMol object (Protein residue)

    returns: The charge of an OBMol object (Protein residue)
    """

    mol.ConnectTheDots()
    mol.PerceiveBondOrders()

    charge = 0

    for atom in ob.OBMolAtomIter(mol):
        if atom.GetAtomicNum() == 7:
            charge += (atom.GetExplicitValence() - 3)

        elif atom.GetAtomicNum() == 8:
            charge += (atom.GetExplicitValence() - 2)

        elif atom.GetAtomicNum() == 16:
            charge += (atom.GetExplicitValence() - 2)

    return charge

def Hydrogenate(carbon, mainAtom, mol):

    """
    Adds 3 hydrogens to a molecule that will be bonded to a carbon

    carbon: The carbon atom where the 3 hydrogens will be added to
    mainAtom: The other atom the carbon is attached to
    mol: The molecule that will be hydrogenated

    returns: THe hydrogenated molecule
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

def GetResidueMolecule(protein, residueNumber):

    """
    Gets the OBMol object corresponding to a residue number on ACE2 protein or Spike protein

    protein: 'ace2' or 'spike'
    residueNumber: The residue Number

    returns: The OBMol object representing the molecule
    """

    resMol = ob.OBMol()

    for res in ob.OBResidueIter(protein):
        if res.GetNum() == residueNumber:
            for atom in ob.OBResidueAtomIter(res):
                resMol.AddAtom(atom)
            break

    return resMol

def GetAtomByID(protein, residueNumber, name):

    """
    Returns the OBAtom object given the AtomID (PDB Columns 13-16),
    the protein name (ace2 or spike), and the residue number
    """

    for res in ob.OBResidueIter(protein):
        if res.GetNum() == residueNumber:
            for atom in ob.OBResidueAtomIter(res):
                if res.GetAtomID(atom) == name:
                    return atom

def NTerminus(protein, residueNumber):
    """
    Returns the N-terminus of a residue given protein name (ace2 or spike),
    and the residue number
    """
    print("executing NTerm")
    return GetAtomByID(protein, residueNumber, ' N  ')

def CTerminus(protein, residueNumber):
    """
    Returns the C-terminus of a residue given protein name (ace2 or spike),
    and the residue number
    """
    print("executing CTerm")
    return GetAtomByID(protein, residueNumber, ' C  ')

def AlphaCarbon(protein, residueNumber):
    """
    Returns the AlphaCarbon of a residue given protein name (ace2 or spike),
    and the residue number
    """
    print("executing AlphaCarbon")
    return GetAtomByID(protein, residueNumber, ' CA ')

def Cap(protein, resNumber):
    """
    Caps the N-terminus and/or C-terminus of an amino acid residue given its protein name (ace2 or spike) and its residue number

    The rules are defined below:
    N-terminus is capped with COCH3 group
    C-terminus is capped with HNCH3 group

    if the residue is at the end of the atom, only the C-terminus and/or N-terminus will be capped, depending on the location of the residue
    """

    resLeft = GetResidueMolecule(protein, resNumber - 1)
    resMol = GetResidueMolecule(protein, resNumber)
    resRight = GetResidueMolecule(protein, resNumber + 1)

    def CapC(resMol):
        resMol.AddAtom(NTerminus(protein, resNumber + 1))
        resMol.AddAtom(AlphaCarbon(protein, resNumber + 1))

        hasH = False

        for atom in ob.OBAtomAtomIter(NTerminus(protein, resNumber + 1)):
            if atom.GetAtomicNum() == 1:
                resMol.AddAtom(atom)
                hasH = True

        if not hasH:
            atom = GetAtomByID(protein, resNumber + 1, ' CD ')

            if atom is not None:

                atom.SetAtomicNum(1)

                origin = np.array([NTerminus(protein, resNumber + 1).GetX(), NTerminus(protein, resNumber + 1).GetY(), NTerminus(protein, resNumber + 1).GetZ()], dtype=float)

                vector = np.array([atom.GetX(), atom.GetY(), atom.GetZ()], dtype=float) - origin

                vector = 1.01*vector/np.linalg.norm(vector)

                vector = origin + vector

                atom.SetVector(vector[0], vector[1], vector[2])

                resMol.AddAtom(atom)

        resMol = Hydrogenate(AlphaCarbon(protein, resNumber + 1), NTerminus(protein, resNumber + 1), resMol)

    def CapN(resMol):
        resMol.AddAtom(CTerminus(protein, resNumber - 1))
        resMol.AddAtom(AlphaCarbon(protein, resNumber - 1))
        resMol.AddAtom(GetAtomByID(protein, resNumber - 1, ' O  '))

        resMol = Hydrogenate(AlphaCarbon(protein, resNumber - 1), CTerminus(protein, resNumber - 1), resMol)


    if resLeft.NumAtoms() == 0:
        CapC(resMol)

    elif resRight.NumAtoms() == 0:
        CapN(resMol)

    else:
        CapN(resMol)
        CapC(resMol)

    return resMol
