import os
from rdkit import Chem
import time
from datetime import timedelta
import json
import paths


start_time = time.time()
input_dir = paths.SEPARATED_SDF_DIR
output_dir = paths.INCHI_DIR
output_json = paths.ALL_INCHIS_DICT
summary_dir = paths.TIMER_DIR
summary_file = os.path.join(summary_dir, "summary_sdf_to_inchi.txt")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(summary_dir):
    os.makedirs(summary_dir)

cofactors = {
    "Ascorbic acid": ["ASC", "UU3"],
    "Coenzyme F420": ["6J4", "F42", "FO1"],
    "Factor F430": ["F43", "M43"],
    "MIO": ["MDO"],
    "Phosphopantetheine": ["PNS"],
    "Pantetheine": ["PNY"],
    "Pantothenic acids": ["66S", "8Q1", "PAU"],
    "Nicotinamide": ["NCA"],
    "Nicotinamide-adenine dinucleotide": [
        "0WD","1DG","3AA","3CD","5J8","6V0","80F","8ID","A3D","AP0",
        "CNA","CND","DG1","DN4","DND","DQV","EAD","ENA","LNC",
        "N01","NA0","NAD","NAE","NAI","NAJ","NAP","NAQ","NAX",
        "NBD","NBP","NDA","NDC","NDE","NDO","NDP","NHD",
        "NHO","NJP","NPW","ODP","P1H","PAD","SAD","SAE","SND",
        "TAD","TAP","TDT","TXD","TXE","TXP","ZID"
    ],
    "Adenosine nucleotides": ["AMP", "ATP", "ADP", "LA8"], # added ADP and LA8 to the list
    "Guanosine nucleotides": ["GTP", "GDP"],
    "Dipyrromethane": ["18W", "29P", "DPM"],
    "Molybdopterin": ["2MD","MCN","MGD","MSS","MTE","MTQ","MTV","PCD","PGD","XAX"],
    "Adenosylcobalamin": ["B12","B1M","CNC","COB","COY"],
    "Flavin adenine dinucleotide": [
        "6FA","FA8","FAA","FAB","FAD","FAE","FAO","FAS","FCG",
        "FDA","FED","FNK","FSH","P5F","RFL","SFD"
    ],
    "Tetrahydrofolic acid": ["1YJ","C2F","DHF","FFO","FOL","FON","FOZ","THF","THG","THH","MEF"],
    "Coenzyme A": [
        "01A","01K","0ET","1C4","1CV","1CZ","1HA","1VU","1XE","2CP","2NE","3CP","3H9","3HC",
        "3VV","4CA","4CO","52O","7L1","8JD","8Z2","94Q","ACO","AMX","BCA","BCO","BSJ","BYC",
        "CA3","CA5","CA6","CA8","CAA","CAJ","CAO","CIC","CMC","CMX","CO6","CO7","CO8","COA",
        "COD","COF","COO","COT","COW","COZ","DCA","DCC","FAM","FCX","FRE","FYN","GRA","HAX",
        "HMG","HSC","HXC","IVC","MCA","MCD","MDE","MLC","MYA","NHM","NHQ","NHW","NMX","OXK",
        "QHD","RMW","S0N","SCA","SCD","SCO","SDX","SOP","T1G","TC6","TUY","WCA","YNC","ZOZ"
    ],
    "Coenzyme B": ["SHT","TP7","TPZ","TXZ","XP8","XP9"],
    "Flavin Mononucleotide": [
        "4LS","4LU","9O9","9P3","9PF","9Q6","9QF","F7F","FMN","FNR","FNS","IRF","RBF"
    ],
    "Lumazine": ["DLZ"],
    "Menaquinone": ["MQ7","MQ8","MQ9","MQE"],
    "Coenzyme M": ["COM"],
    "Heme": [
        "1FH","2FH","522","6HE","76R","7HE","BW9","CCH","COH","CV0","DDH","DHE","F0L","F0X",
        "FDD","FDE","FEC","FMI","H02","HAS","HDD","HDE","HEA","HEB","HEC","HEM","HEO","HEV",
        "HIF","HP5","ISW","MH0","MI9","MNH","MNR","MP1","N7H","OBV","PP9","SH0","SIR","SRM",
        "UFE","VEA","VER","VOV","ZEM","ZNH","1CP","CP3","MMP","UP2","UP3"
    ],
    "Biopterin": ["4AB","7AP","BHS","BIO","H2B","H4B","HBI","WSD"],
    "Methanopterin": ["H4M","H4Z"],
    "Pyrroloquinoline Quinone": ["PQQ"],
    "Biotin": ["BC4","BTI","BTN","BYT","DTB","Y7Y"],
    "Lipoic acid": ["LPA","LPB"],
    "Lipoamide": ["LPM"],
    "Ubiquinone": ["4YP","9BL","AT5","DBT","RQX","UHD","UQ1","UQ2","UQ5","UQ6","DCQ","HQE","PLQ"],
    "Glutathione": [
        "0HG","0HH","1JO","1JP","1R4","3GC","48T","5AU","6SG","ABY","AHE","ATA","BOB","BWS",
        "BYG","EPY","ESG","GBI","GBP","GBX","GDN","GDS","GF5","GGC","GIP","GNB","GPR","GPS",
        "GS8","GSB","GSF","GSH","GSM","GSN","GSO","GTB","GTD","GTS","GTX","GTY","GVX","HAG",
        "HGS","IBG","ICY","JM2","JM5","JM7","L9X","LEE","LZ6","P9H","RGE","TGG","TS5","VWW",
        "ZBF"
    ],
    "Orthoquinone residues (LTQ, TTQ, CTQ)": ["0AF","TOQ","TQQ","TRQ"],
    "S-adenosylmethionine": [
        "0UM","0XU","0Y0","0Y1","0Y2","36A","37H","4IK","62X","6D6","6NR","76H","76J","76K",
        "76L","76M","AN6","EEM","K15","P2J","SA8","SAH","SAM","SFG","SMM","SX0","TT8"
    ],
    "Thiamine diphosphate": [
        "1TP","1U0","2TP","5GY","5SR","8EF","8EL","8EO","8FL","8ML","8N9","8PA","A5X","D7K",
        "EN0","HTL","M6T","N1T","N3T","NDQ","O2T","QSP","R1T","S1T","T5X","T6F","TD5","TD6",
        "TD7","TD8","TD9","TDK","TDL","TDM","TDN","TDP","TDW","THD","THV","THW","THY","TOG",
        "TOI","TP8","TPP","TPU","TPW","TZD","WWF","ZP1"
    ],
    "Pyridoxal": ["PXL","UEG"],
    "Pyridoxal 5'-phosphate": ["EM2","MPL","NOP","NPL","PDP","PLP","PLR","PMP","PXP","PZP","UAH","X04"],
    "Topaquinone": [
        "1TY","2TY","3TY","4HL","AGQ","ESB","G27","HCC","P2Q","P3Q","PAQ","T0I","TPQ","TTS",
        "TYQ","TYY","YPZ"
    ],
    "Siderophores": ["488","EB4","SE8"],
    "Methanofuran": ["MFN"],
    "Vitamin A": ["BCR","ECH","EQ3","RAW"],
    "Vitamin K1": ["PQN"]
}
cofactor_codes = {code for codes in cofactors.values() for code in codes}

def merge_sdf_files(input_dir, output_file):
    # if os.path.exists(output_file):
    #     print(f"Merged SDF file already exists. Skipping merge step.")
    #     return None
    print(f"Merging SDF files ...")
    all_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".sdf")]
    if not all_files:
        print("No SDF files found in input directory.")
        return 0
    writer = Chem.SDWriter(output_file)
    total_mols = 0
    fail_mol = 0
    fail_mol_names = []
    for file_path in all_files:
        sdf_name = os.path.splitext(os.path.basename(file_path))[0]
        suppl = Chem.SDMolSupplier(file_path, removeHs=False)
        for mol in suppl:
            if mol is None:
                fail_mol += 1
                fail_mol_names.append(sdf_name)
                continue
            try:
                mol.SetProp("_Name", sdf_name)
            except Exception:
                pass  
            writer.write(mol)
            total_mols += 1
    writer.close()
    failed_to_merge_count = len(all_files) - total_mols
    # Write the failed molecule names to a txt file in output_dir
    if fail_mol_names:
        failed_txt = os.path.join(output_dir, "failed_to_read_mol_names_from_separated_sdf.txt")
        with open(failed_txt, "w") as fout:
            for name in fail_mol_names:
                fout.write(f"{name}\n")
    return len(all_files), total_mols, failed_to_merge_count, fail_mol

def workflow(input_dir, output_dir):
    # 1. Merge sdf files and Load merged SDF
    merged_sdf = os.path.join(output_dir, "merged_ligands.sdf")
    if os.path.exists(merged_sdf):
        print(f"Merged SDF file already exists. Skipping merge step.")
        total_separated_files, total_mols_merged, failed_to_merge_count, fail_mol = 0, 0, 0, 0
    else:
        print(f"Merged SDF file does not exist. Merging SDF files ...")
        total_separated_files, total_mols_merged, failed_to_merge_count, fail_mol = merge_sdf_files(input_dir, merged_sdf)
    suppl = Chem.SDMolSupplier(merged_sdf, removeHs=False)
    total = 0
    failed = 0
    failed_inchi = 0
    filtered_out_count = 0
    records = []
    for mol in suppl:
        if mol is None:
            failed += 1
            continue
        total += 1
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{total}"
        
        # 2. Apply cofactor filter based on the second token of the name
        identifier = name
        parts = identifier.split("_")
        code = parts[1].upper() if len(parts) >= 2 else ""
        if code in cofactor_codes:
            filtered_out_count += 1
            continue
        inchi_string = Chem.inchi.MolToInchi(mol)
        if inchi_string is None:
            failed_inchi += 1
            continue
        records.append({"Name": name, "InChi": inchi_string})

    # 3. Write JSON dictionary to file
    inchi_dict = {rec["Name"]: rec["InChi"] for rec in records}
    with open(output_json, "w") as f:
        json.dump(inchi_dict, f, indent=4)

    converted = len(records)
    
    return {
        "total_separated_files": total_separated_files,
        "total_mols_merged": total_mols_merged,
        "failed_to_merge_count": failed_to_merge_count,
        "fail_mol": fail_mol,
        "total": total,
        "failed": failed,
        "failed_inchi": failed_inchi,
        "filtered_out": filtered_out_count,
        "converted": converted,
    }

summary_stats = workflow(input_dir, output_dir)

end_time = time.time()
elapsed = end_time - start_time
readable_time = str(timedelta(seconds=int(elapsed)))
print(f"\nTotal script execution time: {readable_time}")

# Persist summary to timer directory
summary_lines = [
    f"Total separated SDF files: {summary_stats.get('total_separated_files', 0)}",
    f"Total molecules merged: {summary_stats.get('total_mols_merged', 0)}",
    f"Failed to merge: {summary_stats.get('failed_to_merge_count', 0)}",
    f"Failed to read molecules from separated SDF files: {summary_stats.get('fail_mol', 0)}",
    f"Total molecules read: {summary_stats.get('total', 0)}",
    f"Failed to parse from merged SDF: {summary_stats.get('failed', 0)}",
    f"Failed to convert to InChI: {summary_stats.get('failed_inchi', 0)}",
    f"Filtered out (cofactors): {summary_stats.get('filtered_out', 0)}",
    f"Converted to InChI: {summary_stats.get('converted', 0)}",
    f"Total script execution time: {readable_time}",
]

for line in summary_lines:
    print(line)

with open(summary_file, "w") as f:
    f.write("\n".join(summary_lines))
print(f"Summary written to {summary_file}")