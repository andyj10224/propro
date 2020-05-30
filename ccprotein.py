from analyze import Cap, GetCharge
import subprocess, os
from openbabel import openbabel as ob
import sys

args = sys.argv

obconv = ob.OBConversion()
obconv.SetInAndOutFormats("pdb", "xyz")
protein = ob.OBMol()

pdb_code = args[1].upper()
chain_letter = args[2].upper()

subprocess.call(['mkdir', F'{pdb_code}'])
os.chdir(F'/Users/andyjiang/Documents/sherrill/CCProtein/{pdb_code}')
subprocess.call(['wget', F'https://files.rcsb.org/download/{pdb_code}.pdb'])

system = open(F"{pdb_code}.pdb")
lines = system.readlines()
system.close()

newfile = open(F'{pdb_code}_{chain_letter}.pdb', 'w+')

for i in range(len(lines)):
    line = lines[i]

    if line[0:6] == "ATOM  " and line[16:17] != 'B':
        if line[21:22] == chain_letter:
            newfile.write(line)

newfile.close()
obconv.ReadFile(protein, F'{pdb_code}_{chain_letter}.pdb')

os.chdir('/Users/andyjiang/Documents/sherrill/CCProtein/')

for res in ob.OBResidueIter(protein):
    mol = Cap(protein, res.GetNum())
    print(res.GetNum(), GetCharge(mol), mol.NumAtoms())
    obconv.WriteFile(mol, F"{pdb_code}/{pdb_code}_{chain_letter}_{res.GetNum()}_capped.xyz")

    resFile = open(F"{pdb_code}/{pdb_code}_{chain_letter}_{res.GetNum()}_capped.xyz", "r")
    resLines = resFile.readlines()
    resFile.close()

    newResFile = open(F"{pdb_code}/{pdb_code}_{chain_letter}_{res.GetNum()}_capped.xyz", "w")
    newResFile.write(resLines[0])
    newResFile.write(F'{GetCharge(mol)} 1\n')

    for i in range(2, len(resLines)):
        newResFile.write(resLines[i])

    newResFile.close()
