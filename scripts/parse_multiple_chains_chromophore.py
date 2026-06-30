import argparse


def main(args):
    import json
    import glob
    import os
    import numpy as np

    folder_with_pdbs_path = args.input_path
    save_path = args.output_path
    ca_only = args.ca_only

    # =========================
    # amino acid mapping
    # =========================
    alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
    states = len(alpha_1)

    alpha_3 = [
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        "GAP",
    ]

    aa_1_N = {a: n for n, a in enumerate(alpha_1)}
    aa_3_N = {a: n for n, a in enumerate(alpha_3)}
    aa_N_1 = {n: a for n, a in enumerate(alpha_1)}
    aa_1_3 = {a: b for a, b in zip(alpha_1, alpha_3)}

    def N_to_AA(x):
        x = np.array(x)
        if x.ndim == 1:
            x = x[None]
        return ["".join([aa_N_1.get(a, "-") for a in y]) for y in x]

    def parse_PDB_biounits(x, atoms=("N", "CA", "C"), chain=None):
        """
        Chromophore-aware parser for fluorescent proteins.

        Compared with original ProteinMPNN parser:
        1. HETATM MSE is converted to MET.
        2. HETATM MLY is converted to LYS.
        3. HETATM CRO / CR2 / NRQ is approximately expanded into
           chromophore-forming parent tripeptide positions.

        Important:
        This does not reconstruct the mature chromophore chemistry.
        It only maps available chromophore atoms to approximate N/CA/C/O
        backbone coordinates for ProteinMPNN input.
        """

        # =========================
        # Chromophore restore rules
        # key: PDB file stem, without .pdb
        #
        # tuple: (chain_id, chromophore_center_resseq, modified_residue_name)
        # parent_seq: restored standard tripeptide sequence
        # offsets: target PDB residue numbers relative to center_resseq
        # =========================
        CHROMO_RULES_BY_PDB = {
            "4EUL": {
                ("A", 66, "CRO"): {
                    "parent_seq": "TYG",
                    "offsets": [-1, 0, 1],
                }
            },
        }

        # Chromophore atom names -> approximate parent tripeptide backbone atoms.
        CHROMO_ATOM_GROUPS = [
            {"N": "N1", "CA": "CA1", "C": "C1", "O": "O1"},
            {"N": "N2", "CA": "CA2", "C": "C2", "O": "O2"},
            {"N": "N3", "CA": "CA3", "C": "C3", "O": "O3"},
        ]

        pdb_name = os.path.basename(x)
        if pdb_name.endswith(".pdb"):
            pdb_name = pdb_name[:-4]

        chromo_rules = (
            CHROMO_RULES_BY_PDB.get(pdb_name, {})
            or CHROMO_RULES_BY_PDB.get(pdb_name.upper(), {})
        )

        xyz = {}
        seq = {}
        min_resn = 1e6
        max_resn = -1e6

        lines = []
        chromo_atoms = {}

        # =========================
        # First pass:
        # collect lines and chromophore HETATM atoms
        # =========================
        with open(x, "rb") as f:
            for raw in f:
                line = raw.decode("utf-8", "ignore").rstrip("\n")
                lines.append(line)

                rec = line[:6]
                resname = line[17:20].strip()
                ch = line[21:22]
                resseq_str = line[22:26].strip()

                if rec.startswith("HETATM") and resname in {"CRO", "CR2", "NRQ"}:
                    if chain is not None and ch != chain:
                        continue

                    try:
                        resseq = int(resseq_str)
                        coord = np.array([float(line[i:(i + 8)]) for i in [30, 38, 46]])
                    except ValueError:
                        continue

                    atom = line[12:16].strip()
                    key = (ch, resseq, resname)
                    chromo_atoms.setdefault(key, {})[atom] = coord

        def add_atom(resn_internal, resa, resi, atom, coord):
            nonlocal min_resn, max_resn

            if resn_internal < min_resn:
                min_resn = resn_internal
            if resn_internal > max_resn:
                max_resn = resn_internal

            if resn_internal not in xyz:
                xyz[resn_internal] = {}
            if resa not in xyz[resn_internal]:
                xyz[resn_internal][resa] = {}

            if resn_internal not in seq:
                seq[resn_internal] = {}
            if resa not in seq[resn_internal]:
                seq[resn_internal][resa] = resi

            if atom not in xyz[resn_internal][resa]:
                xyz[resn_internal][resa][atom] = coord

        # =========================
        # Second pass:
        # parse standard ATOM records and selected modified residues
        # =========================
        for line in lines:
            resname = line[17:20].strip()

            # MSE: selenomethionine -> MET
            if line[:6] == "HETATM" and resname == "MSE":
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace("MSE", "MET")

            # MLY: modified lysine -> LYS
            if line[:6] == "HETATM" and resname == "MLY":
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace("MLY", "LYS")

            if line[:4] != "ATOM":
                continue

            ch = line[21:22]
            if chain is not None and ch != chain:
                continue

            atom = line[12:16].strip()
            resi = line[17:20].strip()
            resn = line[22:27].strip()

            if not resn:
                continue

            try:
                coord = np.array([float(line[i:(i + 8)]) for i in [30, 38, 46]])
            except ValueError:
                continue

            try:
                if resn[-1].isalpha():
                    resa = resn[-1]
                    resn_int = int(resn[:-1]) - 1
                else:
                    resa = ""
                    resn_int = int(resn) - 1
            except ValueError:
                continue

            add_atom(
                resn_internal=resn_int,
                resa=resa,
                resi=resi,
                atom=atom,
                coord=coord,
            )

        # =========================
        # Third step:
        # expand chromophore residues into approximate parent tripeptide backbone
        # =========================
        for (ch, center_resseq, resname), atom_dict in chromo_atoms.items():
            rule = chromo_rules.get((ch, center_resseq, resname))
            if rule is None:
                continue

            parent_seq = rule["parent_seq"]
            offsets = rule.get("offsets", [-1, 0, 1])

            if len(parent_seq) != 3:
                raise ValueError(
                    f"{pdb_name} {ch}{center_resseq} {resname}: parent_seq must have length 3"
                )
            if len(offsets) != 3:
                raise ValueError(
                    f"{pdb_name} {ch}{center_resseq} {resname}: offsets must have length 3"
                )

            print(
                f"[CHROMO] {pdb_name}: {resname} {ch}{center_resseq} "
                f"-> {parent_seq} at offsets {offsets}"
            )

            for i_parent, offset in enumerate(offsets):
                parent_resseq = center_resseq + offset
                resn_internal = parent_resseq - 1

                parent_aa1 = parent_seq[i_parent]
                parent_aa3 = aa_1_3[parent_aa1]

                atom_map = CHROMO_ATOM_GROUPS[i_parent]

                for out_atom in atoms:
                    chromo_atom = atom_map.get(out_atom)
                    if chromo_atom is None:
                        continue
                    if chromo_atom not in atom_dict:
                        continue

                    add_atom(
                        resn_internal=resn_internal,
                        resa="",
                        resi=parent_aa3,
                        atom=out_atom,
                        coord=atom_dict[chromo_atom],
                    )

        # =========================
        # Convert to arrays and fill missing residues/atoms with NaN
        # =========================
        seq_ = []
        xyz_ = []

        try:
            for resn in range(min_resn, max_resn + 1):
                if resn in seq:
                    for k in sorted(seq[resn]):
                        seq_.append(aa_3_N.get(seq[resn][k], 20))
                else:
                    seq_.append(20)

                if resn in xyz:
                    for k in sorted(xyz[resn]):
                        for atom in atoms:
                            if atom in xyz[resn][k]:
                                xyz_.append(xyz[resn][k][atom])
                            else:
                                xyz_.append(np.full(3, np.nan))
                else:
                    for atom in atoms:
                        xyz_.append(np.full(3, np.nan))

            return np.array(xyz_).reshape(-1, len(atoms), 3), N_to_AA(np.array(seq_))

        except TypeError:
            return "no_chain", "no_chain"

    # =========================
    # Main parsing workflow:
    # parse all PDB files and write jsonl
    # =========================
    pdb_dict_list = []

    if folder_with_pdbs_path[-1] != "/":
        folder_with_pdbs_path = folder_with_pdbs_path + "/"

    out_parent = os.path.dirname(os.path.abspath(save_path))
    if out_parent:
        os.makedirs(out_parent, exist_ok=True)

    init_alphabet = [
        "A", "B", "C", "D", "E", "F", "G",
        "H", "I", "J", "K", "L", "M", "N",
        "O", "P", "Q", "R", "S", "T", "U",
        "V", "W", "X", "Y", "Z",
        "a", "b", "c", "d", "e", "f", "g",
        "h", "i", "j", "k", "l", "m", "n",
        "o", "p", "q", "r", "s", "t", "u",
        "v", "w", "x", "y", "z",
    ]
    extra_alphabet = [str(item) for item in list(np.arange(300))]
    chain_alphabet = init_alphabet + extra_alphabet

    biounit_names = sorted(glob.glob(folder_with_pdbs_path + "*.pdb"))

    print("========== parse_multiple_chains_chromophore ==========")
    print(f"[INFO] input_path  = {folder_with_pdbs_path}")
    print(f"[INFO] output_path = {save_path}")
    print(f"[INFO] ca_only     = {ca_only}")
    print(f"[INFO] found PDBs  = {len(biounit_names)}")

    if len(biounit_names) == 0:
        print(f"[WARN] No .pdb files found in: {folder_with_pdbs_path}")

    for biounit in biounit_names:
        my_dict = {}
        s = 0
        concat_seq = ""

        for letter in chain_alphabet:
            if ca_only:
                sidechain_atoms = ["CA"]
            else:
                sidechain_atoms = ["N", "CA", "C", "O"]

            xyz, seq = parse_PDB_biounits(
                biounit,
                atoms=sidechain_atoms,
                chain=letter,
            )

            if type(xyz) != str:
                concat_seq += seq[0]
                my_dict["seq_chain_" + letter] = seq[0]

                coords_dict_chain = {}

                if ca_only:
                    coords_dict_chain["CA_chain_" + letter] = xyz.tolist()
                else:
                    coords_dict_chain["N_chain_" + letter] = xyz[:, 0, :].tolist()
                    coords_dict_chain["CA_chain_" + letter] = xyz[:, 1, :].tolist()
                    coords_dict_chain["C_chain_" + letter] = xyz[:, 2, :].tolist()
                    coords_dict_chain["O_chain_" + letter] = xyz[:, 3, :].tolist()

                my_dict["coords_chain_" + letter] = coords_dict_chain
                s += 1

        fi = biounit.rfind("/")
        my_dict["name"] = biounit[(fi + 1):-4]
        my_dict["num_of_chains"] = s
        my_dict["seq"] = concat_seq

        if s < len(chain_alphabet):
            pdb_dict_list.append(my_dict)

        print(f"[PARSED] {my_dict['name']} | chains={s} | total_length={len(concat_seq)}")

    with open(save_path, "w") as f:
        for entry in pdb_dict_list:
            f.write(json.dumps(entry) + "\n")

    print("========== DONE ==========")
    print(f"[OK] wrote entries = {len(pdb_dict_list)}")
    print(f"[OK] output jsonl  = {save_path}")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument(
        "--input_path",
        type=str,
        help="Path to a folder with pdb files, e.g. /home/my_pdbs/",
    )
    argparser.add_argument(
        "--output_path",
        type=str,
        help="Path where to save .jsonl dictionary of parsed pdbs",
    )
    argparser.add_argument(
        "--ca_only",
        action="store_true",
        default=False,
        help="parse a backbone-only structure (default: false)",
    )

    args = argparser.parse_args()
    main(args)
