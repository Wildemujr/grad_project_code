class BackBone():
    back_bone_dict = {
        "m": ("-\u03B3", "6d"),
        "c": ("-\u03C0", "63"),
        "l": ("-\u03B1", "6c"),
        "p": ("-\u03C1", "70"),
        "b": ("-\u03B2", "62"),
        "t": ("\u03B2", "74"),
        "y": ("\u03C1", "79"),
        "a": ("\u03B1", "61"),
        "i": ("\u03C0", "69"),
        "g": ("\u03B3", "67"),
        "n": ("-\u03B3"+"i", "6e"),
        "d": ("-\u03C0"+"i", "64"),
        "q": ("-\u03B1"+"i", "71"),
        "r": ("-\u03C1"+"i", "72"),
        "f": ("-\u03B2"+"i", "66"),
        "h": ("\u03B2"+"i", "68"),
        "w": ("\u03C1"+"i", "77"),
        "k": ("\u03B1"+"i", "6b"),
        "s": ("\u03C0"+"i", "73"),
        "v": ("\u03B3"+"i", "76"),
        }
    bb_dict_keys = list(back_bone_dict.keys())
    bb_dict_hex_values = [i[1] for i in back_bone_dict.values()]
    bb_dict_greek_values = [i[0] for i in back_bone_dict.values()]
    scoring_matrix = {
        "-γ": {
            "-γ": 7, "-π": 4, "-α": 0, "-ρ": 0, "-β": -1,  "β": 0, "ρ": 0, 
            "α": -1, "π": 3, "γ": 3, "-γi": 4, "-πi": 0, "-αi": 0, "-ρi": 0, 
            "-βi": 0, "βi": 0, "ρi": -2, "αi": 0, "πi": 2, "γi":0
            },
        "-π": {
            "-γ": 4, "-π": 6, "-α": 3, "-ρ": 1, "-β": 0,  "β": 0, "ρ": -1, 
            "α": -1, "π": 1, "γ": 6, "-γi": 8, "-πi": 4, "-αi": 0, "-ρi": -1, 
            "-βi": 3, "βi": 3, "ρi": 0, "αi": 0, "πi": -1, "γi":-1
            },
        "-α": {
            "-γ": 0, "-π": 3, "-α": 4, "-ρ": 2, "-β": 0,  "β": 0, "ρ": -1, 
            "α": -1, "π": 1, "γ": 6, "-γi": 8, "-πi": 4, "-αi": 0, "-ρi": -1, 
            "-βi": 3, "βi": 3, "ρi": 2, "αi": 2, "πi": 2, "γi":1
            },
        "-ρ": {
            "-γ": 0, "-π": 3, "-α": 4, "-ρ": 2, "-β": 0,  "β": 0, "ρ": -1, 
            "α": -1, "π": 1, "γ": 6, "-γi": 8, "-πi": 4, "-αi": 0, "-ρi": -1, 
            "-βi": 3, "βi": 3, "ρi": 2, "αi": 2, "πi": 2, "γi":1
            },
    }


class Alignment(BackBone):
    def hex_to_backbone_alph(self, fasta_seq):
        for index, ele in enumerate(fasta_seq):
            if ele in BackBone.bb_dict_hex_values:
                pos = BackBone.bb_dict_hex_values.index(ele)
                fasta_seq[index] = BackBone.bb_dict_greek_values[pos]
                yield fasta_seq[index]
            elif ele == "--":
                fasta_seq[index] = "---"
                yield fasta_seq[index]





# From protCAD:
# Backbone Classification Type: -γ -π -α -ρ -β β ρ α π γ -γi -πi -αi -ρi -βi βi ρi αi πi γi
# Backbone Classification Key:   m  c  l  p  b t y a i g  n   d   q   r   f  h  w  k  s  v
# Backbone Hex Key: ['6d', '63', '6c', '70', '62', '74', '79', '61', '69', '67', '6e', '64', '71', '72', '66', '68', '77', '6b', '73', '76']
