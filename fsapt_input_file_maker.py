from openbabel import openbabel as ob
#from advanced_charge_finder import CalcTotalCharge
import sys

def main(args=sys.argv):
    template = open("template.in", "r")
    template_lines = template.readlines()
    template.close()

    residue_file = open(args[1], "r")
    residue_lines = residue_file.readlines()
    residue_file.close()

    ligand_file = open(args[2], "r")
    ligand_lines = ligand_file.readlines()
    ligand_file.close()

    conv = ob.OBConversion()

    residue = ob.OBMol()
    ligand = ob.OBMol()
    conv.ReadFile(residue, args[1])
    conv.ReadFile(ligand, args[2])

    residue_charge = args[4]
    ligand_charge = args[5]

    print(residue_charge)
    print(ligand_charge)

    residue_lines[1] = F'{residue_charge} {1}\n'
    ligand_lines[1] = F'{ligand_charge} {1}\n'

    fsapt_input = open(args[3], "w+")

    i = 0

    while (i < len(template_lines)):
        if "RESIDUE" not in template_lines[i] and "LIGAND" not in template_lines[i]:
            fsapt_input.write(template_lines[i])
        else:
            if "RESIDUE_CHARGE" in template_lines[i]:
                fsapt_input.write(residue_lines[1])
            elif "RESIDUE_COORDS" in template_lines[i]:
                for j in range(2, len(residue_lines)):
                    fsapt_input.write(residue_lines[j])
            elif "LIGAND_CHARGE" in template_lines[i]:
                fsapt_input.write(ligand_lines[1])
            elif "LIGAND_COORDS" in template_lines[i]:
                for j in range(2, len(ligand_lines)):
                    fsapt_input.write(ligand_lines[j])

        i += 1

    fsapt_input.close()

if __name__ == '__main__':
    main()
