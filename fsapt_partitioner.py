from openbabel import openbabel as ob
import obaminoacid as obaa
import os

def fsapt_file_creator(res_path, lig_path, fsapt_path):

    obconv = ob.OBConversion()
    residue = ob.OBMol()
    ligand = ob.OBMol()
    obconv.ReadFile(residue, res_path)
    obconv.ReadFile(ligand, lig_path)

    for atom in ob.OBMolAtomIter(residue):
        print(atom.GetIdx(), atom.GetX(), atom.GetY(), atom.GetZ())

    parser = obaa.OBAminoAcidParser(residue)

    NCapIds = []
    CCapIds = []
    Others = []

    for atom in ob.OBMolAtomIter(residue):
        if obaa.OBAminoAcid.IsInMolecule(atom, parser.NTermCap):
            NCapIds.append(atom.GetIdx())

        elif obaa.OBAminoAcid.IsInMolecule(atom, parser.CTermCap):
            CCapIds.append(atom.GetIdx())

        else:
            Others.append(atom.GetIdx())

    fA = open(os.path.join(fsapt_path, 'fA.dat'), 'w+')
    fB = open(os.path.join(fsapt_path, 'fB.dat'), 'w+')

    fA.write("NTermCap ")
    for i in range(len(NCapIds)):
        if i != len(NCapIds) - 1:
            fA.write(F'{NCapIds[i]} ')
        else:
            fA.write(F'{NCapIds[i]}\n')

    fA.write("CTermCap ")
    for i in range(len(CCapIds)):
        if i != len(CCapIds) - 1:
            fA.write(F'{CCapIds[i]} ')
        else:
            fA.write(F'{CCapIds[i]}\n')

    fA.write("Others ")
    for i in range(len(Others)):
        if i != len(Others) - 1:
            fA.write(F'{Others[i]} ')
        else:
            fA.write(F'{Others[i]}')

    fA.close()

    fB.write("Ligand ")

    for atom in ob.OBMolAtomIter(ligand):
        if atom.GetIdx() != ligand.NumAtoms():
            fB.write(F'{residue.NumAtoms() + atom.GetIdx()} ')
        else:
            fB.write(F'{residue.NumAtoms() + atom.GetIdx()}')

    fB.close()
